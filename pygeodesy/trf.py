
# -*- coding: utf-8 -*-

u'''I{Veness}' Terrestrial Reference Frames (TRF).

Classes L{RefFrame}, registry L{RefFrames} and L{TRFError}.

Transcoded from I{Chris Veness'} (C) 2006-2022 JavaScript originals
U{latlon-ellipsoidal-referenceframe.js<https://GitHub.com/ChrisVeness/geodesy/blob/master/
latlon-ellipsoidal-referenceframe.js>} and U{latlon-ellipsoidal-referenceframe-txparams.js
<https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal-referenceframe-txparams.js>}.

Following is a copy of the comments in I{Veness}' U{latlon-ellipsoidal-referenceframe.js
<https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal-referenceframe.js>}.

Modern geodetic reference frames: a latitude/longitude point defines a geographic location on,
above or below the earth’s surface, measured in degrees from the equator and the U{International
Reference Meridian<https://WikiPedia.org/wiki/IERS_Reference_Meridian>} (IRM) and metres above
the ellipsoid within a given I{Terrestrial Reference Frame} at a given I{epoch}.

This is scratching the surface of complexities involved in high precision geodesy, but may
be of interest and/or value to those with less demanding requirements.  More information U{here
<https://www.Movable-Type.co.UK/scripts/geodesy-library.html>} and U{here
<https://www.Movable-Type.co.UK/scripts/geodesy-library.html#latlon-ellipsoidal-referenceframe>}.

Note that I{ITRF solutions} do not directly use an ellipsoid, but are specified by Cartesian
coordinates.  The GRS80 ellipsoid is recommended for transformations to geographical coordinates.

Note WGS84(G730/G873/G1150) are coincident with ITRF at 10-centimetre level, see also U{here
<ftp://ITRF.ENSG.IGN.Fr/pub/itrf/WGS84.TXT>}.  WGS84(G1674) and ITRF20014 / ITRF2008 I{"are likely
to agree at the centimeter level"}, see also U{QPS/Qinsy<https://Confluence.QPS.NL/qinsy/
en/how-to-deal-with-etrs89-datum-and-time-dependent-transformation-parameters-45353274.html>}.

@var RefFrames.ETRF2000: RefFrame(name='ETRF2000', epoch=2005, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.GDA2020: RefFrame(name='GDA2020', epoch=2020, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.GDA94: RefFrame(name='GDA94', epoch=1994, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF2000: RefFrame(name='ITRF2000', epoch=1997, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF2005: RefFrame(name='ITRF2005', epoch=2000, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF2008: RefFrame(name='ITRF2008', epoch=2005, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF2014: RefFrame(name='ITRF2014', epoch=2010, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF90: RefFrame(name='ITRF90', epoch=1988, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF91: RefFrame(name='ITRF91', epoch=1988, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.ITRF93: RefFrame(name='ITRF93', epoch=1988, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.NAD83: RefFrame(name='NAD83', epoch=1997, ellipsoid=Ellipsoid(name='GRS80')
@var RefFrames.WGS84: RefFrame(name='WGS84', epoch=1984, ellipsoid=Ellipsoid(name='WGS84')
@var RefFrames.WGS84g1150: RefFrame(name='WGS84g1150', epoch=2001, ellipsoid=Ellipsoid(name='WGS84')
@var RefFrames.WGS84g1674: RefFrame(name='WGS84g1674', epoch=2005, ellipsoid=Ellipsoid(name='WGS84')
@var RefFrames.WGS84g1762: RefFrame(name='WGS84g1762', epoch=2005, ellipsoid=Ellipsoid(name='WGS84')
'''

from pygeodesy.basics import map1, _xinstanceof, _zip
from pygeodesy.constants import _float as _F, _0_0s, _0_0, _0_001, \
                                _0_01, _0_1, _0_26, _0_5, _N_0_5, _1_0
from pygeodesy.datums import _a_ellipsoid, Transform
from pygeodesy.ellipsoids import Ellipsoids, Fmt, Property_RO
from pygeodesy.errors import _ALL_LAZY, _IsnotError, TRFError
from pygeodesy.interns import NN, _COMMASPACE_, _cartesian_, _conversion_, \
                             _ellipsoid_, _ellipsoidal_, _epoch_, _exists_, \
                             _GRS80_, _NAD83_, _name_, _no_, _s_, _SPACE_, \
                             _sx_, _sy_, _sz_, _to_, _tx_, _ty_, _tz_, _WGS84_
# from pygeodesy.lazily import _ALL_LAZY  # from .errors
from pygeodesy.named import classname, _lazyNamedEnumItem as _lazy, _NamedEnum, \
                           _NamedEnumItem, _NamedDict as _XD, _NamedTuple
# from pygeodesy.props import Property_RO  # from .ellipsoids
# from pygeodesy.streprs import Fmt  # from .ellipsoids
from pygeodesy.units import Epoch, Float

from math import ceil

__all__ = _ALL_LAZY.trf
__version__ = '23.05.26'

_0_02  =   0.02
_0_06  =   0.06
_0_09  =   0.09
_0_12  =   0.12
_366_0 = 366.0

_Forward =  _0_001  # mm2m, ppb2ppM, mas2as
_Inverse = -_0_001  # same, inverse transforms
_mas     =  _mm = _ppb = Float  # as == arcseconds
_mDays   =  (0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0)
_N_0_06  = -_0_06
_N_0_1   = -_0_1

_ETRF2000_   = 'ETRF2000'
_GDA2020_    = 'GDA2020'
_GDA94_      = 'GDA94'
_ITRF_       = 'ITRF'
_ITRF88_     = 'ITRF88'
_ITRF89_     = 'ITRF89'
_ITRF90_     = 'ITRF90'
_ITRF91_     = 'ITRF91'
_ITRF92_     = 'ITRF92'
_ITRF93_     = 'ITRF93'
_ITRF94_     = 'ITRF94'
_ITRF96_     = 'ITRF96'
_ITRF97_     = 'ITRF97'
_ITRF2000_   = 'ITRF2000'
_ITRF2005_   = 'ITRF2005'
_ITRF2008_   = 'ITRF2008'
_ITRF2014_   = 'ITRF2014'
_WGS84g1150_ = 'WGS84g1150'
_WGS84g1674_ = 'WGS84g1674'
_WGS84g1762_ = 'WGS84g1762'


class Helmert7Tuple(_NamedTuple):
    '''7-Tuple C{(tx, ty, tz, s, sx, sy, sz)} of Helmert transform
       parameters with translations C{tx}, C{ty} and C{tz} in
       C{millimeter}, scale C{s} in C{ppb} and rotations C{sx},
       C{sy} and C{sz} in C{milliarcseconds}.

       @see: L{Datum} L{Transform}.
    '''
    _Names_ = (_tx_, _ty_, _tz_, _s_,  _sx_, _sy_, _sz_)
    _Units_ = (_mm,  _mm,  _mm,  _ppb, _mas, _mas, _mas)

    def __new__(cls, tx, ty, tz, s, sx, sy, sz, name=NN):
        '''New L{Helmert7Tuple}.
        '''
        t = map1(_F, tx, ty, tz, s, sx, sy, sz)
        return _NamedTuple.__new__(cls, *t, name=name)


class RefFrame(_NamedEnumItem):
    '''Terrestrial Reference Frame (TRF) parameters.
    '''
    _ellipsoid =  None  # ellipsoid GRS80 or WGS84 (L{Ellipsoid} or L{Ellipsoid2})
    _epoch     = _0_0   # epoch, calendar year (L{Epoch} or C{float})

    def __init__(self, epoch, ellipsoid, name=NN):
        '''New L{RefFrame}.

           @arg epoch: Epoch, a fractional calendar year (C{scalar} or C{str}).
           @arg ellipsoid: The ellipsoid (L{Ellipsoid}, L{Ellipsoid2},
                           L{Datum} or L{a_f2Tuple}).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: A L{RefFrame} with that B{C{name}}
                             already exists.

           @raise TRFError: Invalid B{C{epoch}}.

           @raise TypeError: Invalid B{C{ellipsoid}}.
        '''
        self._ellipsoid = _a_ellipsoid(ellipsoid, name=name, raiser=_ellipsoid_)
        self._epoch = Epoch(epoch)
        self._register(RefFrames, name)

    def __matmul__(self, other):  # PYCHOK Python 3.5+
        '''Convert cartesian or ellipsoidal B{C{other}} to this reframe.

           @raise TypeError: Invalid B{C{other}}.
        '''
        try:  # only Cartesian- and LatLonEllipsoidalBase
            return other.toRefFrame(self)
        except AttributeError:
            pass
        raise _IsnotError(_cartesian_, _ellipsoidal_, other=other)

    @Property_RO
    def ellipsoid(self):
        '''Get this reference frame's ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return self._ellipsoid

    @Property_RO
    def epoch(self):
        '''Get this reference frame's epoch (C{Epoch}).
        '''
        return self._epoch

    def toStr(self, name=NN, **unused):  # PYCHOK expected
        '''Return this reference frame as a text string.

           @kwarg name: Override name (C{str}) or C{None} to exclude the
                        reframe's name.

           @return: This L{RefFrame}'s attributes (C{str}).
        '''
        e = self.ellipsoid
        t = (Fmt.EQUAL(_name_, repr(name or self.name)),
             Fmt.EQUAL(_epoch_, self.epoch),
             Fmt.PAREN(Fmt.EQUAL(_ellipsoid_, classname(e)),
                       Fmt.EQUAL(_name_, repr(e.name))))
        return _COMMASPACE_.join(t[1:] if name is None else t)


class RefFrames(_NamedEnum):
    '''(INTERNAL) L{RefFrame} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, epoch, ellipsoid_name, name=NN):
        '''(INTERNAL) Instantiate the L{RefFrame}.
        '''
        return RefFrame(epoch, Ellipsoids.get(ellipsoid_name), name=name)

RefFrames = RefFrames(RefFrame)  # PYCHOK singleton
'''Some pre-defined L{RefFrame}s, all I{lazily} instantiated.'''
# <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal-referenceframe.js>
RefFrames._assert(
    ETRF2000   = _lazy(_ETRF2000_,   _F(2005), _GRS80_),  # ETRF2000(R08)
    GDA2020    = _lazy(_GDA2020_,    _F(2020), _GRS80_),  # Australia
    GDA94      = _lazy(_GDA94_,      _F(1994), _GRS80_),  # Australia
    ITRF2000   = _lazy(_ITRF2000_,   _F(1997), _GRS80_),
    ITRF2005   = _lazy(_ITRF2005_,   _F(2000), _GRS80_),
    ITRF2008   = _lazy(_ITRF2008_,   _F(2005), _GRS80_),  # aka ITRF08
    ITRF2014   = _lazy(_ITRF2014_,   _F(2010), _GRS80_),
    ITRF90     = _lazy(_ITRF90_,     _F(1988), _GRS80_),
    ITRF91     = _lazy(_ITRF91_,     _F(1988), _GRS80_),
    ITRF93     = _lazy(_ITRF93_,     _F(1988), _GRS80_),
    NAD83      = _lazy(_NAD83_,      _F(1997), _GRS80_),  # CORS96
    WGS84      = _lazy(_WGS84_,      _F(1984), _WGS84_),
    WGS84g1150 = _lazy(_WGS84g1150_, _F(2001), _WGS84_),
    WGS84g1674 = _lazy(_WGS84g1674_, _F(2005), _WGS84_),
    WGS84g1762 = _lazy(_WGS84g1762_, _F(2005), _WGS84_))  # same epoch


def date2epoch(year, month, day):
    '''Return the reference frame C{epoch} for a calendar day.

       @arg year: Year of the date (C{scalar}).
       @arg month: Month in the B{C{year}} (C{scalar}, 1..12).
       @arg day: Day in the B{C{month}} (C{scalar}, 1..31).

       @return: Epoch, the fractional year (C{float}).

       @raise TRFError: Invalid B{C{year}}, B{C{month}} or B{C{day}}.

       @note: Any B{C{year}} is considered a leap year, i.e. having
              29 days in February.
    '''
    try:
        y, m, d = map1(int, year, month, day)
        if y > 0 and 1 <= m <= 12 and 1 <= d <= _mDays[m]:
            return Epoch(y + float(sum(_mDays[:m]) + d) / _366_0, low=0)

        raise ValueError  # _invalid_
    except (TRFError, TypeError, ValueError) as x:
        raise TRFError(year=year, month=month, day=day, cause=x)


def epoch2date(epoch):
    '''Return the date for a reference frame C{epoch}.

       @arg epoch: Fractional year (C{scalar}).

       @return: 3-Tuple C{(year, month, day)}.

       @raise TRFError: Invalid B{C{epoch}}.

       @note: Any B{C{year}} is considered a leap year, i.e. having
              29 days in February.
    '''
    e = Epoch(epoch, Error=TRFError, low=0)
    y = int(e)
    d = int(ceil(_366_0 * (e - y)))
    for m, n in enumerate(_mDays[1:]):
        if d > n:
            d -= n
        else:
            break
    return y, (m + 1), max(1, d)


def _intermediate(n1, n2):
    '''(INTERNAL) Find a trf* "in between" C{n1} and C{n2}.
    '''
    f1 = set(m for n, m in _trfXs.keys() if n == n1)  # from trf1
    t2 = set(n for n, m in _trfXs.keys() if m == n2)  # to trf2
    n  = f1.intersection(t2)
    return n.pop() if n else NN


def _reframeTransforms2(rf2, rf, epoch):
    '''(INTERNAL) Get 0, 1 or 2 Helmert L{Transform}s to convert
       reference frame B{C{rf}} observed at B{C{epoch}} into B{C{rf2}}.
    '''
    e = rf.epoch if epoch is None else Epoch(epoch)

    n2 = rf2.name  # .upper()
    n1 = rf.name   # .upper()
    if n1 == n2 or (n1.startswith(_ITRF_) and n2.startswith(_WGS84_)) \
                or (n2.startswith(_ITRF_) and n1.startswith(_WGS84_)):
        return e, ()  # PYCHOK returns

    if (n1, n2) in _trfXs:
        return e, (_2Transform((n1, n2), e, _Forward),)  # PYCHOK returns

    if (n2, n1) in _trfXs:
        return e, (_2Transform((n2, n1), e, _Inverse),)  # PYCHOK returns

    n = _intermediate(n1, n2)
    if n:
        return e, (_2Transform((n1, n), e, _Forward),  # PYCHOK returns
                   _2Transform((n, n2), e, _Forward))

    n = _intermediate(n2, n1)
    if n:
        return e, (_2Transform((n, n1), e, _Inverse),  # PYCHOK returns
                   _2Transform((n2, n), e, _Inverse))

    t = _SPACE_(RefFrame.__name__, repr(n1), _to_, repr(n2))
    raise TRFError(_no_(_conversion_), txt=t)


def _2Transform(n1_n2, epoch, _Forward_Inverse):
    '''(INTERNAL) Combine two Helmert transforms from TRF
       conversion C{_trfXs[n1_n2]} observed at B{C{epoch}}
       into a single datum L{Transform}.

       @note: Translations are converted from C{millimeter}
              to C{meter}, rotations from C{milliarcseconds}
              to C{arcseconds} and scale from C{ppb} to C{ppM}.
    '''
    X = _trfXs[n1_n2]
    e =  epoch - X.epoch  # fractional delta years
    d =  dict((n, (p + r * e) * _Forward_Inverse) for n, p, r in
              _zip(Helmert7Tuple._Names_, X.xform, X.rates))  # strict=True
    return Transform(**d)


def trfXform(reframe1, reframe2, epoch=None, xform=None, rates=None):
    '''Define a new Terrestrial Reference Frame (TRF) conversion.

       @arg reframe1: Source reframe (L{RefFrame}), converting I{from}.
       @arg reframe2: Destination reframe (L{RefFrame}), converting I{to}.
       @kwarg epoch: Epoch, a fractional calendar year (C{scalar} or C{str})
                     or C{None} for C{B{reframe2}.epoch}.
       @kwarg xform: Helmert transform parameters (C{Helmert7Tuple}).
       @kwarg rates: Helmert rate parameters (C{Helmert7Tuple}).

       @raise TRFError: Invalid B{C{epoch}} or TRF conversion already exists.
    '''
    _xinstanceof(RefFrame, reframe1=reframe1, reframe2=reframe2)
    e = reframe2.epoch if epoch is None else Epoch(epoch=epoch, Error=TRFError)
    _xinstanceof(Helmert7Tuple, xform=xform, rates=rates)
    _trfX(reframe1.name, reframe2.name, epoch=e, xform=xform, rates=rates)


def _trfX(n1, n2, **epoch_xform_rates):
    '''(INTERNAL) New C{_trfXs} entry.
    '''
    n1_n2 = n1, n2
    if n1_n2 in _trfXs:
        raise TRFError(trfX=n1_n2, txt=_exists_)  # _NameError
    _trfXs[n1_n2] = _XD(X=n1_n2, **epoch_xform_rates)


_H      =  Helmert7Tuple
_H_0_0s = _H(*_0_0s(len(_H._Names_)))
# TRF conversions specified as an epoch and dual 7-parameter Helmert transforms.  Most from U{Transformation
# Parameters<http://ITRF.IGN.Fr/trans_para.php>}, more at U{Quinsy QPS<https://confluence.QPS.NL/qinsy/
# files/latest/en/182618383/182618384/1/1579182881000/ITRF_Transformation_Parameters.xlsx>}.  See also
# U{Quinsy International Terrestrial Reference Frame 2014 (ITRF2014)<https://confluence.QPS.NL/qinsy/
# latest/en/international-terrestrial-reference-frame-2014-itrf2014-182618383.html>}.
_trfXs = dict()  # key is [(from_TRF, to_TRF)] 2-tuple
# see U{Transformation Parameters ITRF2014<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>}
_trfX(_ITRF2014_, _ITRF2008_, epoch=_F(2010),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2014/tp_14-08.php>
                              xform=_H(   1.6,      1.9,     2.4,   -0.02,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,  _N_0_1,    0.03,     _0_0,     _0_0,      _0_0))
_trfX(_ITRF2014_, _ITRF2005_, epoch=_F(2010),
                              xform=_H(   2.6,     _1_0,    -2.3,    0.92,     _0_0,     _0_0,      _0_0),
                              rates=_H(   0.3,     _0_0,  _N_0_1,    0.03,     _0_0,     _0_0,      _0_0))
_trfX(_ITRF2014_, _ITRF2000_, epoch=_F(2010),
                              xform=_H(   0.7,      1.2,   -26.1,    2.12,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_1,     _0_1,    -1.9,    0.11,     _0_0,     _0_0,      _0_0))
_trfX(_ITRF2014_, _ITRF97_,   epoch=_F(2010),
                              xform=_H(   7.4,   _N_0_5,   -62.8,    3.8,      _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF96_,   epoch=_F(2010),
                              xform=_H(   7.4,   _N_0_5,   -62.8,    3.8,      _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF94_,   epoch=_F(2010),
                              xform=_H(   7.4,   _N_0_5,   -62.8,    3.8,      _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF93_,   epoch=_F(2010),
                              xform=_H( -50.4,      3.3,   -60.2,    4.29,     -2.81,    -3.38,      0.4),
                              rates=_H(  -2.8,   _N_0_1,    -2.5,   _0_12,     -0.11,    -0.19,      0.07))
_trfX(_ITRF2014_, _ITRF92_,   epoch=_F(2010),
                              xform=_H(  15.4,      1.5,   -70.8,    3.09,     _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF91_,   epoch=_F(2010),
                              xform=_H(  27.4,     15.5,   -76.8,    4.49,     _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF90_,   epoch=_F(2010),
                              xform=_H(  25.4,     11.5,   -92.8,    4.79,     _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF89_,   epoch=_F(2010),
                              xform=_H(  30.4,     35.5,  -130.8,    8.19,     _0_0,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2014_, _ITRF88_,   epoch=_F(2010),
                              xform=_H(  25.4,   _N_0_5,  -154.8,   11.29,     _0_1,     _0_0,      _0_26),
                              rates=_H(  _0_1,   _N_0_5,    -3.3,   _0_12,     _0_0,     _0_0,      _0_02))

# see U{Transformation Parameters ITRF2008<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt>}
# _trfX(_ITRF2008_, _ITRF2005_, epoch=_F(2005),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2008/tp_08-05.php>
#                             xform=_H(_N_0_5,  s  -0.9,    -4.7,    0.94,     _0_0,     _0_0,      _0_0),
#                             rates=_H(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0))
_trfX(_ITRF2008_, _ITRF2005_, epoch=_F(2000),
                              xform=_H(  -2.0,     -0.9,    -4.7,    0.94,     _0_0,     _0_0,      _0_0),
                              rates=_H(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0))
_trfX(_ITRF2008_, _ITRF2000_, epoch=_F(2000),
                              xform=_H(  -1.9,     -1.7,   -10.5,    1.34,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_1,     _0_1,    -1.8,    0.08,     _0_0,     _0_0,      _0_0))
_trfX(_ITRF2008_, _ITRF97_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF96_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF94_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF93_,   epoch=_F(2000),
                              xform=_H( -24.0,      2.4,   -38.6,    3.41,     -1.71,    -1.48,     -0.3),
                              rates=_H(  -2.8,   _N_0_1,    -2.4,   _0_09,     -0.11,    -0.19,      0.07))
_trfX(_ITRF2008_, _ITRF92_,   epoch=_F(2000),
                              xform=_H(  12.8,      4.6,   -41.2,    2.21,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF91_,   epoch=_F(2000),
                              xform=_H(  24.8,     18.6,   -47.2,    3.61,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF90_,   epoch=_F(2000),
                              xform=_H(  22.8,     14.6,   -63.2,    3.91,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF89_,   epoch=_F(2000),
                              xform=_H(  27.8,     38.6,  -101.2,    7.31,     _0_0,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2008_, _ITRF88_,   epoch=_F(2000),
                              xform=_H(  22.8,      2.6,  -125.2,   10.41,     _0_1,     _0_0,      _0_06),
                              rates=_H(  _0_1,   _N_0_5,    -3.2,   _0_09,     _0_0,     _0_0,      _0_02))

_trfX(_ITRF2005_, _ITRF2000_, epoch=_F(2000),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2005/tp_05-00.php>
                              xform=_H(  _0_1,     -0.8,    -5.8,    0.4,      _0_0,     _0_0,      _0_0),
                              rates=_H(  -0.2,     _0_1,    -1.8,    0.08,     _0_0,     _0_0,      _0_0))

_trfX(_ITRF2000_, _ITRF97_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF96_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF94_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF93_,   epoch=_F(1988),
                              xform=_H(  12.7,      6.5,   -20.9,    1.95,     -0.39,     0.8,      -1.14),
                              rates=_H(  -2.9,     -0.2,    -0.6,   _0_01,     -0.11,    -0.19,      0.07))
_trfX(_ITRF2000_, _ITRF92_,   epoch=_F(1988),
                              xform=_H(   1.47,     1.35,   -1.39,   0.75,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF91_,   epoch=_F(1988),
                              xform=_H(   26.7,    27.5,   -19.9,    2.15,     _0_0,     _0_0,      -0.18),
                              rates=_H(   _0_0,    -0.6,    -1.4,   _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF90_,   epoch=_F(1988),
                              xform=_H(   2.47,     2.35,   -3.59,   2.45,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF89_,   epoch=_F(1988),
                              xform=_H(   2.97,     4.75,   -7.39,   5.85,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))
_trfX(_ITRF2000_, _ITRF88_,   epoch=_F(1988),
                              xform=_H(   2.47,     1.15,   -9.79,   8.95,     _0_1,     _0_0,      -0.18),
                              rates=_H(  _0_0,   _N_0_06,   -0.14,  _0_01,     _0_0,     _0_0,      _0_02))

# see U{Boucher, C. & Altamimi, Z. "Memo: Specifications for reference frame fixing in the
# analysis of a EUREF GPS campaign" (2011) <https://ETRS89.ENSG.IGN.Fr/memo-V8.pdf>} and
# Altamimi, Z. U{"Key results of ITRF2014 and implication to ETRS89 realization", EUREF2016
# <https://www.EUREF.EU/symposia/2016SanSebastian/01-02-Altamimi.pdf>}.
_trfX(_ITRF2014_, _ETRF2000_, epoch=_F(2000),
                              xform=_H(  53.7,     51.2,   -55.1,    1.02,      0.891,    5.39,     -8.712),
                              rates=_H(  _0_1,     _0_1,    -1.9,    0.11,      0.081,    0.49,     -0.792))
_trfX(_ITRF2008_, _ETRF2000_, epoch=_F(2000),
                              xform=_H(  52.1,     49.3,   -58.5,    1.34,      0.891,    5.39,     -8.712),
                              rates=_H(  _0_1,     _0_1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2005_, _ETRF2000_, epoch=_F(2000),
                              xform=_H(  54.1,     50.2,   -53.8,    0.4,       0.891,    5.39,     -8.712),
                              rates=_H(  -0.2,     _0_1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2000_, _ETRF2000_, epoch=_F(2000),
                              xform=_H(  54.0,     51.0,   -48.0,   _0_0,       0.891,    5.39,     -8.712),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.792))

# GDA2020 "Geocentric Datum of Australia 2020 Technical Manual", v1.5, 2020-12-09, Table 3.3 and 3.4
# <https://www.ICSM.gov.AU/sites/default/files/2020-12/GDA2020%20Technical%20Manual%20V1.5_4.pdf>
# (the GDA2020 xforms are different but the rates are the same as GDA94, further below)
_trfX(_ITRF2014_, _GDA2020_,  epoch=_F(2020),
                              xform=_H_0_0s,
                              rates=_H(  _0_0,     _0_0,   _0_0,    _0_0,       1.50379,  1.18346,   1.20716))
_trfX(_ITRF2008_, _GDA2020_,  epoch=_F(2020),
                              xform=_H(  13.79,     4.55,   15.22,   2.5,       0.2808,   0.2677,   -0.4638),
                              rates=_H(   1.42,     1.34,    0.9,    0.109,     1.5461,   1.182,     1.1551))
_trfX(_ITRF2005_, _GDA2020_,  epoch=_F(2020),
                              xform=_H(  40.32,   -33.85,  -16.72,   4.286,    -1.2893,  -0.8492,   -0.3342),
                              rates=_H(   2.25,    -0.62,   -0.56,   0.294,    -1.4707,  -1.1443,   -1.1701))
_trfX(_ITRF2000_, _GDA2020_,  epoch=_F(2020),
                              xform=_H(-105.52,    51.58,  231.68,   3.55,      4.2175,   6.3941,    0.8617),
                              rates=_H(  -4.66,     3.55,   11.24,   0.249,     1.7454,   1.4868,    1.224))

# see Table 2 in U{Dawson, J. & Woods, A. "ITRF to GDA94 coordinate transformations", Journal of Applied
# Geodesy 4 (2010), 189-199<https://www.ResearchGate.net/publication/258401581_ITRF_to_GDA94_coordinate_transformations>}
# (note, sign of rotations for GDA94 reversed as "Australia assumes rotation to be of coordinate axes",
# rather than the more conventional "position around the coordinate axes")
_trfX(_ITRF2008_, _GDA94_,    epoch=_F(1994),
                              xform=_H( -84.68,   -19.42,   32.01,   9.71,     -0.4254,   2.2578,    2.4015),
                              rates=_H(   1.42,     1.34,    0.9,    0.109,     1.5461,   1.182,     1.1551))
_trfX(_ITRF2005_, _GDA94_,    epoch=_F(1994),
                              xform=_H( -79.73,    -6.86,   38.03,   6.636,     0.0351,  -2.1211,   -2.1411),
                              rates=_H(   2.25,    -0.62,   -0.56,   0.294,    -1.4707,  -1.1443,   -1.1701))
_trfX(_ITRF2000_, _GDA94_,    epoch=_F(1994),
                              xform=_H( -45.91,   -29.85,  -20.37,   7.07,     -1.6705,   0.4594,    1.9356),
                              rates=_H(  -4.66,     3.55,   11.24,   0.249,     1.7454,   1.4868,    1.224))

# see U{Solar, T. & Snay, R.A. "Transforming Positions and Velocities between the
# International Terrestrial Reference Frame of 2000 and North American Datum of 1983"
# (2004)<https://www.NGS.NOAA.gov/CORS/Articles/SolerSnayASCE.pdf>}
_trfX(_ITRF2000_, _NAD83_,    epoch=_F(1997),  # note NAD83(CORS96)
                              xform=_H( 995.6,  -1901.3,  -521.5,    0.615,   -25.915,   -9.426,   -11.599),
                              rates=_H(   0.7,     -0.7,    _0_5,   -0.182,    -0.06667,  0.75744,   0.05133))

# see U{Quinsy QPS<https://confluence.QPS.NL/qinsy/files/latest/en/182618383/182618384/1/1579182881000/
# ITRF_Transformation_Parameters.xlsx>}, sheets ITRF and NAD83
_trfX(_ITRF2008_, _NAD83_,    epoch=_F(1997),
                              xform=_H( 993.43, -1903.31, -526.55,   1.71504, -25.91467, -9.42645, -11.59935),
                              rates=_H(   0.79,    -0.6,    -1.34,  -0.10201,  -0.06667,  0.75744,   0.05133))
_trfX(_ITRF2005_, _NAD83_,    epoch=_F(1997),
                              xform=_H( 996.3,  -1902.4,  -521.9,    0.775,   -25.915,   -9.426,   -11.599),
                              rates=_H(   0.5,     -0.6,    -1.3,   -0.10201,  -0.06667,  0.75744,   0.05133))
_trfX(_ITRF90_,   _NAD83_,    epoch=_F(1997),
                              xform=_H( 973.0,  -1919.2, -482.9,    -0.9,     -25.79,    -9.65,    -11.66),
                              rates=_H(  _0_0,     _0_0,   _0_0,    _0_0,      -0.053,    0.742,     0.032))
_trfX(_ITRF90_,   _WGS84_,    epoch=_F(1984),
                              xform=_H(  60.0,   -517.0, -223.0,   -11.0,      18.3,     -0.3,       7.0),
                              rates=_H_0_0s)
# del _H, _H_0_0s

if __name__ == '__main__':

    from pygeodesy.interns import _COMMA_, _NL_, _NLATvar_, _STAR_
    from pygeodesy.lazily import printf
    from time import localtime

    D = date2epoch.__name__
    E = epoch2date.__name__
    y = localtime()[0]
    for m in range(1, 13):
        for d in (1, 15, _mDays[m] - 1, _mDays[m]):
            f = '%s(%d,%3d,%3d)' % (D, y, m, d)
            e = date2epoch(y, m, d)
            t = epoch2date(e)
            x = NN if t == (y, m, d) else _STAR_
            e = '%.3f' % (e,)
            e = '%s, %s(%s)' % (e, E, e)
            t = '%d,%3d,%3d' % t
            printf('# %s = %s = %s %s', f, e, t, x)

    # __doc__ of this file, force all into registery
    t = [NN] + RefFrames.toRepr(all=True).split(_NL_)
    printf(_NLATvar_.join(i.strip(_COMMA_) for i in t))

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

# % python -m pygeodesy.trf
#
# date2epoch(2022,  1,  1) = 2022.003, epoch2date(2022.003) = 2022,  1,  1
# date2epoch(2022,  1, 15) = 2022.041, epoch2date(2022.041) = 2022,  1, 15
# date2epoch(2022,  1, 30) = 2022.082, epoch2date(2022.082) = 2022,  1, 30
# date2epoch(2022,  1, 31) = 2022.085, epoch2date(2022.085) = 2022,  1, 31
# date2epoch(2022,  2,  1) = 2022.087, epoch2date(2022.087) = 2022,  2,  2 *
# date2epoch(2022,  2, 15) = 2022.126, epoch2date(2022.126) = 2022,  2, 16 *
# date2epoch(2022,  2, 28) = 2022.161, epoch2date(2022.161) = 2022,  2, 28
# date2epoch(2022,  2, 29) = 2022.164, epoch2date(2022.164) = 2022,  3,  1 *
# date2epoch(2022,  3,  1) = 2022.167, epoch2date(2022.167) = 2022,  3,  2 *
# date2epoch(2022,  3, 15) = 2022.205, epoch2date(2022.205) = 2022,  3, 16 *
# date2epoch(2022,  3, 30) = 2022.246, epoch2date(2022.246) = 2022,  3, 31 *
# date2epoch(2022,  3, 31) = 2022.249, epoch2date(2022.249) = 2022,  4,  1 *
# date2epoch(2022,  4,  1) = 2022.251, epoch2date(2022.251) = 2022,  4,  1
# date2epoch(2022,  4, 15) = 2022.290, epoch2date(2022.290) = 2022,  4, 15
# date2epoch(2022,  4, 29) = 2022.328, epoch2date(2022.328) = 2022,  4, 29
# date2epoch(2022,  4, 30) = 2022.331, epoch2date(2022.331) = 2022,  4, 30
# date2epoch(2022,  5,  1) = 2022.333, epoch2date(2022.333) = 2022,  5,  1
# date2epoch(2022,  5, 15) = 2022.372, epoch2date(2022.372) = 2022,  5, 15
# date2epoch(2022,  5, 30) = 2022.413, epoch2date(2022.413) = 2022,  5, 30
# date2epoch(2022,  5, 31) = 2022.415, epoch2date(2022.415) = 2022,  6,  1 *
# date2epoch(2022,  6,  1) = 2022.418, epoch2date(2022.418) = 2022,  6,  2 *
# date2epoch(2022,  6, 15) = 2022.456, epoch2date(2022.456) = 2022,  6, 16 *
# date2epoch(2022,  6, 29) = 2022.495, epoch2date(2022.495) = 2022,  6, 30 *
# date2epoch(2022,  6, 30) = 2022.497, epoch2date(2022.497) = 2022,  7,  1 *
# date2epoch(2022,  7,  1) = 2022.500, epoch2date(2022.500) = 2022,  7,  1
# date2epoch(2022,  7, 15) = 2022.538, epoch2date(2022.538) = 2022,  7, 16 *
# date2epoch(2022,  7, 30) = 2022.579, epoch2date(2022.579) = 2022,  7, 30
# date2epoch(2022,  7, 31) = 2022.582, epoch2date(2022.582) = 2022,  7, 31
# date2epoch(2022,  8,  1) = 2022.585, epoch2date(2022.585) = 2022,  8,  1
# date2epoch(2022,  8, 15) = 2022.623, epoch2date(2022.623) = 2022,  8, 15
# date2epoch(2022,  8, 30) = 2022.664, epoch2date(2022.664) = 2022,  8, 31 *
# date2epoch(2022,  8, 31) = 2022.667, epoch2date(2022.667) = 2022,  9,  1 *
# date2epoch(2022,  9,  1) = 2022.669, epoch2date(2022.669) = 2022,  9,  2 *
# date2epoch(2022,  9, 15) = 2022.708, epoch2date(2022.708) = 2022,  9, 16 *
# date2epoch(2022,  9, 29) = 2022.746, epoch2date(2022.746) = 2022,  9, 30 *
# date2epoch(2022,  9, 30) = 2022.749, epoch2date(2022.749) = 2022, 10,  1 *
# date2epoch(2022, 10,  1) = 2022.751, epoch2date(2022.751) = 2022, 10,  1
# date2epoch(2022, 10, 15) = 2022.790, epoch2date(2022.790) = 2022, 10, 15
# date2epoch(2022, 10, 30) = 2022.831, epoch2date(2022.831) = 2022, 10, 30
# date2epoch(2022, 10, 31) = 2022.833, epoch2date(2022.833) = 2022, 10, 31
# date2epoch(2022, 11,  1) = 2022.836, epoch2date(2022.836) = 2022, 11,  1
# date2epoch(2022, 11, 15) = 2022.874, epoch2date(2022.874) = 2022, 11, 15
# date2epoch(2022, 11, 29) = 2022.913, epoch2date(2022.913) = 2022, 11, 29
# date2epoch(2022, 11, 30) = 2022.915, epoch2date(2022.915) = 2022, 12,  1 *
# date2epoch(2022, 12,  1) = 2022.918, epoch2date(2022.918) = 2022, 12,  2 *
# date2epoch(2022, 12, 15) = 2022.956, epoch2date(2022.956) = 2022, 12, 16 *
# date2epoch(2022, 12, 30) = 2022.997, epoch2date(2022.997) = 2022, 12, 31 *
# date2epoch(2022, 12, 31) = 2023.000, epoch2date(2023.000) = 2023,  1,  1 *
