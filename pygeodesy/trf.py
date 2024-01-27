
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

@var RefFrames.ETRF88: RefFrame(name='ETRF88', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ETRF89: RefFrame(name='ETRF89', epoch=1989, datum=Datums.GRS80)
@var RefFrames.ETRF90: RefFrame(name='ETRF90', epoch=1990, datum=Datums.GRS80)
@var RefFrames.ETRF91: RefFrame(name='ETRF91', epoch=1991, datum=Datums.GRS80)
@var RefFrames.ETRF92: RefFrame(name='ETRF92', epoch=1992, datum=Datums.GRS80)
@var RefFrames.ETRF93: RefFrame(name='ETRF93', epoch=1993, datum=Datums.GRS80)
@var RefFrames.ETRF94: RefFrame(name='ETRF94', epoch=1994, datum=Datums.GRS80)
@var RefFrames.ETRF96: RefFrame(name='ETRF96', epoch=1996, datum=Datums.GRS80)
@var RefFrames.ETRF97: RefFrame(name='ETRF97', epoch=1997, datum=Datums.GRS80)
@var RefFrames.ETRF2000: RefFrame(name='ETRF2000', epoch=2005, datum=Datums.GRS80)
@var RefFrames.ETRF2005: RefFrame(name='ETRF2005', epoch=2005, datum=Datums.GRS80)
@var RefFrames.ETRF2008: RefFrame(name='ETRF2008', epoch=2008, datum=Datums.GRS80)
@var RefFrames.ETRF2014: RefFrame(name='ETRF2014', epoch=2014, datum=Datums.GRS80)
@var RefFrames.GDA94: RefFrame(name='GDA94', epoch=1994, datum=Datums.GRS80)
@var RefFrames.GDA2020: RefFrame(name='GDA2020', epoch=2020, datum=Datums.GRS80)
@var RefFrames.ITRF88: RefFrame(name='ITRF88', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ITRF89: RefFrame(name='ITRF89', epoch=1989, datum=Datums.GRS80)
@var RefFrames.ITRF90: RefFrame(name='ITRF90', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ITRF91: RefFrame(name='ITRF91', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ITRF92: RefFrame(name='ITRF92', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ITRF93: RefFrame(name='ITRF93', epoch=1988, datum=Datums.GRS80)
@var RefFrames.ITRF94: RefFrame(name='ITRF94', epoch=1993, datum=Datums.GRS80)
@var RefFrames.ITRF96: RefFrame(name='ITRF96', epoch=1997, datum=Datums.GRS80)
@var RefFrames.ITRF97: RefFrame(name='ITRF97', epoch=1997, datum=Datums.GRS80)
@var RefFrames.ITRF2000: RefFrame(name='ITRF2000', epoch=1997, datum=Datums.GRS80)
@var RefFrames.ITRF2005: RefFrame(name='ITRF2005', epoch=2000, datum=Datums.GRS80)
@var RefFrames.ITRF2008: RefFrame(name='ITRF2008', epoch=2005, datum=Datums.GRS80)
@var RefFrames.ITRF2014: RefFrame(name='ITRF2014', epoch=2010, datum=Datums.GRS80)
@var RefFrames.ITRF2020: RefFrame(name='ITRF2020', epoch=2015, datum=Datums.GRS80)
@var RefFrames.NAD83: RefFrame(name='NAD83', epoch=1997, datum=Datums.GRS80)
@var RefFrames.WGS84: RefFrame(name='WGS84', epoch=1984, datum=Datums.GRS80)
@var RefFrames.WGS84g1150: RefFrame(name='WGS84g1150', epoch=2001, datum=Datums.GRS80)
@var RefFrames.WGS84g1674: RefFrame(name='WGS84g1674', epoch=2005, datum=Datums.GRS80)
@var RefFrames.WGS84g1762: RefFrame(name='WGS84g1762', epoch=2005, datum=Datums.GRS80)
'''

from pygeodesy.basics import map1, isstr, _xinstanceof, _zip
from pygeodesy.constants import _float as _F, _0_0, _0_001, _0_5, _1_0
from pygeodesy.datums import Datums, _earth_datum, Transform, _WGS84
from pygeodesy.errors import _IsnotError, TRFError
from pygeodesy.interns import MISSING, NN, _AT_, _COMMASPACE_, _cartesian_, _conversion_, \
                             _datum_, _DOT_, _ellipsoidal_, _exists_, _NAD83_, _no_, \
                             _reframe_, _s_, _SPACE_, _sx_, _sy_, _sz_, _to_, _tx_, _ty_, \
                             _tz_, _WGS84_, _intern as _i
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import classname, _lazyNamedEnumItem as _lazy, _NamedDict as _XD, \
                           _NamedEnum, _NamedEnumItem, _NamedTuple,  Fmt
from pygeodesy.props import Property_RO, property_RO
# from pygeodesy.streprs import Fmt  # from .named
from pygeodesy.units import Epoch, Float

from math import ceil

__all__ = _ALL_LAZY.trf
__version__ = '24.01.25'

_366_0   = 366.0
_Forward =  _0_001  # mm2m, ppb2ppM, mas2as
_GRS80   =   Datums.GRS80
_Inverse = -_0_001  # same, inverse transforms
_mas     =  _mm = _ppb = Float  # as == arcseconds
_mDays   =  (0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0)

_ETRF88_     = _i('ETRF88')
_ETRF89_     = _i('ETRF89')
_ETRF90_     = _i('ETRF90')
_ETRF91_     = _i('ETRF91')
_ETRF92_     = _i('ETRF92')
_ETRF93_     = _i('ETRF93')
_ETRF94_     = _i('ETRF94')
_ETRF96_     = _i('ETRF96')
_ETRF97_     = _i('ETRF97')
_ETRF2000_   = _i('ETRF2000')
_ETRF2005_   = _i('ETRF2005')
_ETRF2008_   = _i('ETRF2008')
_ETRF2014_   = _i('ETRF2014')
_GDA94_      = _i('GDA94')
_GDA2020_    = _i('GDA2020')
_ITRF_       = _i('ITRF')
_ITRF88_     = _i('ITRF88')
_ITRF89_     = _i('ITRF89')
_ITRF90_     = _i('ITRF90')
_ITRF91_     = _i('ITRF91')
_ITRF92_     = _i('ITRF92')
_ITRF93_     = _i('ITRF93')
_ITRF94_     = _i('ITRF94')
_ITRF96_     = _i('ITRF96')
_ITRF97_     = _i('ITRF97')
_ITRF2000_   = _i('ITRF2000')
_ITRF2005_   = _i('ITRF2005')
_ITRF2008_   = _i('ITRF2008')
_ITRF2014_   = _i('ITRF2014')
_ITRF2020_   = _i('ITRF2020')
_WGS84g1150_ = _i('WGS84g1150')
_WGS84g1674_ = _i('WGS84g1674')
_WGS84g1762_ = _i('WGS84g1762')
# del _i


class Helmert7Tuple(_NamedTuple):
    '''7-Tuple C{(tx, ty, tz, s, sx, sy, sz)} of Helmert transform
       parameters with translations C{tx}, C{ty} and C{tz} in
       C{millimeter}, scale C{s} in C{ppb} and rotations C{sx},
       C{sy} and C{sz} in C{milli-arc-seconds}.

       @note: The parameters names are often capitalized and
              alternate names are C{(Tx, Ty, Tz, D, Rx, Ry, Rz)}
              and C{(T1, T2, T3, D, R1, R2, R3)}.

       @see: L{Datum} L{Transform} keyword arguments.
    '''
    _Names_ = (_tx_, _ty_, _tz_, _s_,  _sx_, _sy_, _sz_)  # == kwds
    _Units_ = (_mm,  _mm,  _mm,  _ppb, _mas, _mas, _mas)

    def __new__(cls, tx, ty, tz, s, sx, sy, sz, name=NN):
        '''New L{Helmert7Tuple}.
        '''
        t = map1(_F, tx, ty, tz, s, sx, sy, sz)
        return _NamedTuple.__new__(cls, *t, name=name)


class RefFrame(_NamedEnumItem):
    '''Terrestrial Reference Frame (TRF) parameters.
    '''
    _datum = _GRS80  # Datums.GRS80 or .WGS84 (L{Datum})
    _epoch = _0_0    # epoch, calendar year (L{Epoch} or C{float})

    def __init__(self, epoch, datum, name=NN):
        '''New L{RefFrame}.

           @arg epoch: Epoch, a fractional calendar year (C{scalar} or C{str}).
           @arg datum: Datum or ellipsoid (L{Datum}, {Ellipsoid}, L{Ellipsoid2}
                       or L{a_f2Tuple}).
           @kwarg name: Unique, non-empty name (C{str}).

           @raise NameError: A L{RefFrame} with that B{C{name}} already exists.

           @raise TRFError: Invalid B{C{epoch}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        if datum is not _GRS80:
            _earth_datum(self, datum, raiser=_datum_)
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

    @property_RO
    def datum(self):
        '''Get this reference frame's datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def ellipsoid(self):
        '''Get this reference frame's ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return self._datum.ellipsoid

    @Property_RO
    def epoch(self):
        '''Get this reference frame's epoch (C{Epoch}).
        '''
        return self._epoch

    def toRefFrame(self, point, reframe2, epoch2=None, epoch=None, name=NN):
        '''Convert a cartesian or geodetic point from this to another reframe and epoch.

           @return: A copy of the B{C{point}}, converted or renamed.

           @see: Ellipsoidal methods L{LatLon.toRefFrame<ellipsoidalBase.LatLonEllipsoidalBase.toRefFrame>}
                 and L{Cartesian.toRefFrame<ellipsoidalBase.CartesianEllipsoidalBase.toRefFrame>}
                 for more details.
        '''
        b = _MODS.ellipsoidalBase
        _xinstanceof(b.LatLonEllipsoidalBase, b.CartesianEllipsoidalBase, point=point)
        r = point.dup(reframe=self)
        return r.toRefFrame(reframe2, epoch2=epoch2, epoch=epoch, name=name or self.name)

    def toStr(self, epoch=None, name=NN, **unused):  # PYCHOK expected
        '''Return this reference frame as a text string.

           @kwarg epoch: Override this reframe's epoch (C{scalar} or C{str}).
           @kwarg name: Override name (C{str}) or C{None} to exclude the
                        reframe's name.

           @return: This L{RefFrame}'s attributes (C{str}).
        '''
        D = self.datum
        e = self.epoch if epoch is None else Epoch(epoch)
        t = (Fmt.EQUAL(name=repr(name or self.name)),
             Fmt.EQUAL(epoch=e),
             Fmt.EQUAL(datum=_DOT_(classname(D) + _s_, D.name)))
        return _COMMASPACE_.join(t[1:] if name is None else t)


class RefFrames(_NamedEnum):
    '''(INTERNAL) L{RefFrame} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, epoch, datum=_GRS80, name=NN):
        '''(INTERNAL) Instantiate the L{RefFrame}.
        '''
        return RefFrame(epoch, datum, name=name)

RefFrames = RefFrames(RefFrame)  # PYCHOK singleton
'''Some pre-defined L{RefFrame}s, all I{lazily} instantiated.'''
# <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal-referenceframe.js>
RefFrames._assert(
    ETRF88     = _lazy(_ETRF88_,     _F(1988)),  # epoch?
    ETRF89     = _lazy(_ETRF89_,     _F(1989)),  # epoch?
    ETRF90     = _lazy(_ETRF90_,     _F(1990)),  # epoch?
    ETRF91     = _lazy(_ETRF91_,     _F(1991)),  # epoch?
    ETRF92     = _lazy(_ETRF92_,     _F(1992)),  # epoch?
    ETRF93     = _lazy(_ETRF93_,     _F(1993)),  # epoch?
    ETRF94     = _lazy(_ETRF94_,     _F(1994)),  # epoch?
    ETRF96     = _lazy(_ETRF96_,     _F(1996)),  # epoch?
    ETRF97     = _lazy(_ETRF97_,     _F(1997)),  # epoch?
    ETRF2000   = _lazy(_ETRF2000_,   _F(2005)),
    ETRF2005   = _lazy(_ETRF2005_,   _F(2005)),  # epoch?
    ETRF2008   = _lazy(_ETRF2008_,   _F(2008)),  # epoch?
    ETRF2014   = _lazy(_ETRF2014_,   _F(2014)),  # epoch?
    GDA94      = _lazy(_GDA94_,      _F(1994)),  # Australia
    GDA2020    = _lazy(_GDA2020_,    _F(2020)),  # Australia
    ITRF88     = _lazy(_ITRF88_,     _F(1988)),
    ITRF89     = _lazy(_ITRF89_,     _F(1989)),
    ITRF90     = _lazy(_ITRF90_,     _F(1988)),
    ITRF91     = _lazy(_ITRF91_,     _F(1988)),
    ITRF92     = _lazy(_ITRF92_,     _F(1988)),
    ITRF93     = _lazy(_ITRF93_,     _F(1988)),
    ITRF94     = _lazy(_ITRF94_,     _F(1993)),
    ITRF96     = _lazy(_ITRF96_,     _F(1997)),
    ITRF97     = _lazy(_ITRF97_,     _F(1997)),
    ITRF2000   = _lazy(_ITRF2000_,   _F(1997)),
    ITRF2005   = _lazy(_ITRF2005_,   _F(2000)),
    ITRF2008   = _lazy(_ITRF2008_,   _F(2005)),  # aka ITRF08
    ITRF2014   = _lazy(_ITRF2014_,   _F(2010)),
    ITRF2020   = _lazy(_ITRF2020_,   _F(2015)),
    NAD83      = _lazy(_NAD83_,      _F(1997)),  # CORS96
    WGS84      = _lazy(_WGS84_,      _F(1984), _WGS84),
    WGS84g1150 = _lazy(_WGS84g1150_, _F(2001), _WGS84),
    WGS84g1674 = _lazy(_WGS84g1674_, _F(2005), _WGS84),
    WGS84g1762 = _lazy(_WGS84g1762_, _F(2005), _WGS84))  # same epoch


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
    e = Epoch(epoch, low=0)
    y = int(e)
    d = int(ceil(_366_0 * (e - y)))
    for m, n in enumerate(_mDays[1:]):
        if d > n:
            d -= n
        else:
            break
    return y, (m + 1), max(1, d)


def _eTsDs4(inst, reframe, epoch, reframe2, epoch2):
    '''(INTERNAL) Get epoch, a 0-, 1- or 2-tuple of Helmert L{Transform}s
       datum and datum2 to convert B{C{refFrame}} observed at B{C{epoch}}
       into B{C{refFrame2}} observed at B{C{epoch2 or epoch}}.
    '''
    r = reframe or inst.reframe
    if not r:
        t = _SPACE_(_DOT_(repr(inst), _reframe_), MISSING)
        raise TRFError(_no_(_conversion_), txt=t)

    _xinstanceof(RefFrame, reframe2=reframe2, reframe=r)

    e1 = Epoch(epoch or inst.epoch or r.epoch)
    e2 = e1 if epoch2 is None else Epoch(epoch2=epoch2)

    xs = _2Transforms(r.name, e1, reframe2.name, e2)
    if xs is None:
        t = _SPACE_(RefFrame.__name__, repr(r.name), _AT_, e1,
                          _to_, repr(reframe2.name), _AT_, e2)
        raise TRFError(_no_(_conversion_), txt=t)

    return e2, xs, r.datum, reframe2.datum


def _intermediate(n1, n2):
    '''(INTERNAL) Find a trf* "between" C{n1} and C{n2}.
    '''
    f = set(m for n, m in _trfXs.keys() if n == n1)  # from trf1
    t = set(n for n, m in _trfXs.keys() if m == n2)  # to trf2
    n = f.intersection(t)
    return n.pop() if n else NN


def _2Transform(n1, n2, e, _Forward_Inverse):
    '''(INTERNAL) Combine the TRF C{xform} and C{rates} from a
       conversion C{_trfXs[(n1, n2)]} observed at C{B{e}poch}
       into a single I{datum} L{Transform}.

       @note: Translations are converted from C{millimeter} to C{meter},
              rotations from C{milliarcseconds} to C{arcseconds} and
              scale from C{ppb} to C{ppM}.
    '''
    X = _trfXs[(n1, n2)]
    e -= X.epoch  # delta in fractional years
    d  = dict((n, (x + e * r) * _Forward_Inverse) for n, x, r in
              _zip(Helmert7Tuple._Names_, X.xform, X.rates))  # strict=True
    return Transform(**d)


def _2Transforms(n1, e1, n2, e2):
    '''(INTERNAL) Get 0-, 1- or 2-tuple of Helmert L{Transform}s or C{None}.
    '''
    if n1 == n2 or (n1.startswith(_ITRF_) and n2.startswith(_WGS84_)) \
                or (n2.startswith(_ITRF_) and n1.startswith(_WGS84_)):
        return ()  # PYCHOK returns
    if (n1, n2) in _trfXs:
        return (_2Transform(n1, n2, e1, _Forward),  # PYCHOK returns
                _2Transform(n1, n2, e2, _Forward)) if e2 != e1 else \
               (_2Transform(n1, n2, e2, _Forward),)
    if (n2, n1) in _trfXs:
        return (_2Transform(n2, n1, e1, _Inverse),  # PYCHOK returns
                _2Transform(n2, n1, e2, _Inverse)) if e2 != e1 else \
               (_2Transform(n2, n1, e2, _Inverse),)
    n = _intermediate(n1, n2)
    if n:
        return (_2Transform(n1, n, e1, _Forward),  # PYCHOK returns
                _2Transform(n, n2, e2, _Forward))
    n = _intermediate(n2, n1)
    if n:
        return (_2Transform(n, n1, e1, _Inverse),  # PYCHOK returns
                _2Transform(n2, n, e2, _Inverse))
    return None


def trfTransforms(reframe, epoch, reframe2, epoch2):
    '''Get the L{Transform}(s) to convert reframe at epoch to reframe2 at epoch2.

       @arg reframe: Reference frame to convert I{from} (L{RefFrame} or C{str}).
       @arg epoch: Epoch to observe I{from} (L{Epoch}, C{scalar} or C{str}).
       @arg reframe2: Reference frame to convert I{to} (L{RefFrame} or C{str}).
       @arg epoch2: Epoch to observe to observe I{to} (L{Epoch}, C{scalar} or C{str}).

       @return: A 0-, 1- or 2-tuple of Helmert L{Transform}s or C{None} if no
                conversion exists.
    '''
    _xinstanceof(RefFrame, str, reframe=reframe, reframe2=reframe2)
    n1 = reframe.upper()  if isstr(reframe)  else reframe.name
    n2 = reframe2.upper() if isstr(reframe2) else reframe2.name
    return _2Transforms(n1, Epoch(epoch), n2, Epoch(epoch2=epoch2))


def trfXform(reframe1, reframe2, epoch=None, xform=None, rates=None):
    '''Define a new Terrestrial Reference Frame (TRF) conversion.

       @arg reframe1: Source reframe (L{RefFrame}), converting I{from}.
       @arg reframe2: Destination reframe (L{RefFrame}), converting I{to}.
       @kwarg epoch: Epoch, a fractional calendar year (C{scalar} or C{str})
                     or C{None} for C{B{reframe2}.epoch}.
       @kwarg xform: Helmert I{transform} parameters (C{Helmert7Tuple}).
       @kwarg rates: Helmert I{rate} parameters (C{Helmert7Tuple}), like
                     B{C{xform}}, but in C{units per year}.

       @raise TRFError: Invalid B{C{epoch}} or TRF conversion already exists.
    '''
    _xinstanceof(RefFrame, reframe1=reframe1, reframe2=reframe2)
    e = reframe2.epoch if epoch is None else Epoch(epoch=epoch)
    _xinstanceof(Helmert7Tuple, xform=xform, rates=rates)
    _trfX(reframe1.name, reframe2.name, epoch=e, xform=xform, rates=rates)


def _trfX(n1, n2, **epoch_xform_rates):
    '''(INTERNAL) New C{_trfXs} entry.
    '''
    n1_n2 = n1, n2
    if n1_n2 in _trfXs:
        raise TRFError(trfX=n1_n2, txt=_exists_)  # _NameError
    _trfXs[n1_n2] = _XD(X=n1_n2, **epoch_xform_rates)


def _H(*ps):
    h = Helmert7Tuple(*ps)
    return _Hs.setdefault(h, h)  # PYCHOK del

# TRF conversions specified as an epoch and dual 7-parameter Helmert transforms.  Most from U{Transformation
# Parameters<http://ITRF.IGN.Fr/trans_para.php>}, more at U{Quinsy QPS<https://confluence.QPS.NL/qinsy/
# files/latest/en/182618383/182618384/1/1579182881000/ITRF_Transformation_Parameters.xlsx>}.  See also
# U{Quinsy International Terrestrial Reference Frame 2014 (ITRF2014)<https://confluence.QPS.NL/qinsy/
# latest/en/international-terrestrial-reference-frame-2014-itrf2014-182618383.html>}.
_Hs    = dict()  # PYCHOK unique HelmertTuples, temporary
_trfXs = dict()  # key is [(from_TRF, to_TRF)] 2-tuple
_trfX(_ITRF2020_, _ITRF2014_, epoch=_F(2015),  # <https://ITRF.IGN.Fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt>
                              xform=_H(  -1.4,     -0.9,     1.4,   -0.42,     _0_0,      _0_0,     _0_0),
                              rates=_H(  _0_0,     -0.1,     0.2,    _0_0,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2008_, epoch=_F(2015),
                              xform=_H(   0.2,      1.0,     3.3,   -0.29,     _0_0,      _0_0,     _0_0),
                              rates=_H(  _0_0,     -0.1,     0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2005_, epoch=_F(2015),
                              xform=_H(   2.7,      0.1,    -1.4,    0.65,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.3,     -0.1,     0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2000_, epoch=_F(2015),
                              xform=_H(  -0.2,      0.8,   -34.2,    2.25,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.1,     _0_0,    -1.7,    0.11,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF97_,   epoch=_F(2015),
                              xform=_H(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF96_,   epoch=_F(2015),
                              xform=_H(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF94_,   epoch=_F(2015),
                              xform=_H(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF93_,   epoch=_F(2015),
                              xform=_H( -65.8,      1.9,   -71.3,    4.47,     -3.36,     -4.33,     0.75),
                              rates=_H(  -2.8,     -0.2,    -2.3,    0.12,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2020_, _ITRF92_,   epoch=_F(2015),
                              xform=_H(  14.5,     -1.9,   -85.9,    3.27,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF91_,   epoch=_F(2015),
                              xform=_H(  26.5,     12.1,   -91.9,    4.67,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF90_,   epoch=_F(2015),
                              xform=_H(  24.5,      8.1,  -107.9,    4.97,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF89_,   epoch=_F(2015),
                              xform=_H(  29.5,     32.1,  -145.9,    8.37,     _0_0,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF88_,   epoch=_F(2015),
                              xform=_H(  24.5,     -3.9,  -169.9,   11.47,      0.1,      _0_0,      0.36),
                              rates=_H(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))

# see U{Transformation Parameters ITRF2014<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>} and
# Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Appendix A.
_trfX(_ITRF2014_, _ITRF2008_, epoch=_F(2010),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2014/tp_14-08.php>
                              xform=_H(   1.6,      1.9,     2.4,   -0.02,     _0_0,      _0_0,     _0_0),
                              rates=_H(  _0_0,     _0_0,    -0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF2005_, epoch=_F(2010),
                              xform=_H(   2.6,     _1_0,    -2.3,    0.92,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.3,     _0_0,    -0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF2000_, epoch=_F(2010),
                              xform=_H(   0.7,      1.2,   -26.1,    2.12,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.1,      0.1,    -1.9,    0.11,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF97_,   epoch=_F(2010),
                              xform=_H(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF96_,   epoch=_F(2010),
                              xform=_H(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF94_,   epoch=_F(2010),
                              xform=_H(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF93_,   epoch=_F(2010),
                              xform=_H( -50.4,      3.3,   -60.2,    4.29,     -2.81,     -3.38,     0.4),
                              rates=_H(  -2.8,     -0.1,    -2.5,    0.12,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2014_, _ITRF92_,   epoch=_F(2010),
                              xform=_H(  15.4,      1.5,   -70.8,    3.09,     _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF91_,   epoch=_F(2010),
                              xform=_H(  27.4,     15.5,   -76.8,    4.49,     _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF90_,   epoch=_F(2010),
                              xform=_H(  25.4,     11.5,   -92.8,    4.79,     _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF89_,   epoch=_F(2010),
                              xform=_H(  30.4,     35.5,  -130.8,    8.19,     _0_0,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF88_,   epoch=_F(2010),
                              xform=_H(  25.4,     -0.5,  -154.8,   11.29,      0.1,      _0_0,      0.26),
                              rates=_H(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))

# see U{Transformation Parameters ITRF2008<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt>}
# _trfX(_ITRF2008_, _ITRF2005_, epoch=_F(2005),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2008/tp_08-05.php>
#                             xform=_H(  -0.5,  s  -0.9,    -4.7,    0.94,     _0_0,      _0_0,     _0_0),
#                             rates=_H(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,     _0_0))
_trfX(_ITRF2008_, _ITRF2005_, epoch=_F(2000),
                              xform=_H(  -2.0,     -0.9,    -4.7,    0.94,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,     _0_0))
_trfX(_ITRF2008_, _ITRF2000_, epoch=_F(2000),
                              xform=_H(  -1.9,     -1.7,   -10.5,    1.34,     _0_0,      _0_0,     _0_0),
                              rates=_H(   0.1,      0.1,    -1.8,    0.08,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2008_, _ITRF97_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF96_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF94_,   epoch=_F(2000),
                              xform=_H(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF93_,   epoch=_F(2000),
                              xform=_H( -24.0,      2.4,   -38.6,    3.41,     -1.71,     -1.48,    -0.3),
                              rates=_H(  -2.8,     -0.1,    -2.4,    0.09,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2008_, _ITRF92_,   epoch=_F(2000),
                              xform=_H(  12.8,      4.6,   -41.2,    2.21,     _0_0,      _0_0,      0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF91_,   epoch=_F(2000),
                              xform=_H(  24.8,     18.6,   -47.2,    3.61,     _0_0,     _0_0,       0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF90_,   epoch=_F(2000),
                              xform=_H(  22.8,     14.6,   -63.2,    3.91,     _0_0,     _0_0,       0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF89_,   epoch=_F(2000),
                              xform=_H(  27.8,     38.6,  -101.2,    7.31,     _0_0,     _0_0,       0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF88_,   epoch=_F(2000),
                              xform=_H(  22.8,      2.6,  -125.2,   10.41,      0.1,     _0_0,       0.06),
                              rates=_H(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))

_trfX(_ITRF2005_, _ITRF2000_, epoch=_F(2000),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2005/tp_05-00.php>
                              xform=_H(   0.1,     -0.8,    -5.8,    0.4,      _0_0,     _0_0,      _0_0),
                              rates=_H(  -0.2,      0.1,    -1.8,    0.08,     _0_0,     _0_0,      _0_0))

_trfX(_ITRF2000_, _ITRF97_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF96_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF94_,   epoch=_F(1997),
                              xform=_H(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF93_,   epoch=_F(1988),
                              xform=_H(  12.7,      6.5,   -20.9,    1.95,     -0.39,     0.8,      -1.14),
                              rates=_H(  -2.9,     -0.2,    -0.6,    0.01,     -0.11,    -0.19,      0.07))
_trfX(_ITRF2000_, _ITRF92_,   epoch=_F(1988),
                              xform=_H(   1.47,     1.35,   -1.39,   0.75,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF91_,   epoch=_F(1988),
                              xform=_H(   26.7,    27.5,   -19.9,    2.15,     _0_0,     _0_0,      -0.18),
                              rates=_H(   _0_0,    -0.6,    -1.4,    0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF90_,   epoch=_F(1988),
                              xform=_H(   2.47,     2.35,   -3.59,   2.45,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF89_,   epoch=_F(1988),
                              xform=_H(   2.97,     4.75,   -7.39,   5.85,     _0_0,     _0_0,      -0.18),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF88_,   epoch=_F(1988),
                              xform=_H(   2.47,     1.15,   -9.79,   8.95,      0.1,     _0_0,      -0.18),
                              rates=_H(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))

# see Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Table 1.
# _trfX(_ITRF2014_, _ETRF2014_, epoch=_F(1989),
#                             xform=_H(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0),
#                             rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.085,    0.531,    -0.77))
_trfX(_ITRF2005_, _ETRF2005_, epoch=_F(1989),
                              xform=_H(  56.0,     48.0,   -37.0,   _0_0,      _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.054,    0.518,    -0.781))
# _trfX(_ITRF2000_, _ETRF2000_, epoch=_F(1989),
#                             xform=_H(  54.0,     51.0,   -48.0,   _0_0,      _0_0,     _0_0,      _0_0),
#                             rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.792))
_trfX(_ITRF97_,   _ETRF97_,   epoch=_F(1989),
                              xform=_H(  41.0,     41.0,   -49.0,   _0_0,      _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,      0.5,      -0.65))
_trfX(_ITRF96_,   _ETRF96_,   epoch=_F(1989),
                              xform=_H(  41.0,     41.0,   -49.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,       0.5,      -0.65))
_trfX(_ITRF94_,   _ETRF94_,   epoch=_F(1989),
                              xform=_H(  41.0,     41.0,   -49.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,       0.5,      -0.65))
_trfX(_ITRF93_,   _ETRF93_,   epoch=_F(1989),
                              xform=_H(  19.0,     53.0,   -21.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.32,      0.78,     -0.67))
_trfX(_ITRF92_,   _ETRF92_,   epoch=_F(1989),
                              xform=_H(  38.0,     40.0,   -37.0,    0.0,       0.0,       0.0,       0.0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.21,      0.52,     -0.68))
_trfX(_ITRF91_,   _ETRF91_,   epoch=_F(1989),
                              xform=_H(  21.0,     25.0,   -37.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.21,      0.52,     -0.68))
_trfX(_ITRF90_,   _ETRF90_,   epoch=_F(1989),
                              xform=_H(  19.0,     28.0,   -23.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.11,      0.57,     -0.71))
_trfX(_ITRF89_,   _ETRF89_,   epoch=_F(1989),
                              xform=_H(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.11,      0.57,     -0.71))

# see Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Table 2.
_trfX(_ITRF2014_, _ETRF2014_, epoch=_F(2010),
                              xform=_H(  _0_0,     _0_0,    _0_0,   _0_0,       1.785,    11.151,   -16.17),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.085,     0.531,    -0.77))
_trfX(_ITRF2008_, _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -1.6,     -1.9,    -2.4,    0.02,      1.785,    11.151,   -16.17),
                              rates=_H(  _0_0,     _0_0,     0.1,   -0.03,      0.085,     0.531,    -0.77))
_trfX(_ITRF2005_, _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -2.6,     -1.0,     2.3,   -0.92,      1.785,    11.151,   -16.17),
                              rates=_H(  -0.3,     _0_0,     0.1,   -0.03,      0.085,     0.531,    -0.77))
_trfX(_ITRF2000_, _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -0.7,     -1.2,    26.1,   -2.12,      1.785,    11.151,   -16.17),
                              rates=_H(  -0.1,     -0.1,     1.9,   -0.11,      0.085,     0.531,    -0.77))
_trfX(_ITRF97_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -7.4,      0.5,    62.8,   -3.8,       1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF96_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -7.4,      0.5,    62.8,   -3.8,       1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,  -0.12,       0.085,     0.531,    -0.79))
_trfX(_ITRF94_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H(  -7.4,      0.5,    62.8,   -3.8,       1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF93_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H(  50.4,     -3.3,    60.2,   -4.29,      4.595,    14.531,   -16.57),
                              rates=_H(   2.8,      0.1,     2.5,   -0.12,      0.195,     0.721,    -0.84))
_trfX(_ITRF92_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H( -15.4,     -1.5,    70.8,   -3.09,      1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF91_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H( -27.4,    -15.5,    76.8,   -4.49,      1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF90_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H( -25.4,    -11.5,    92.8,   -4.79,      1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF89_,   _ETRF2014_, epoch=_F(2010),
                              xform=_H( -30.4,    -35.5,   130.8,   -8.19,      1.785,    11.151,   -16.43),
                              rates=_H(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))

# see U{Altamimi, Z. "EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Table 3,
# U{Boucher, C. & Altamimi, Z. "Memo: Specifications for reference frame fixing in the analysis of a EUREF GPS
# campaign" (2011) <https://ETRS89.ENSG.IGN.Fr/memo-V8.pdf>} and U{Altamimi, Z. "Key results of ITRF2014 and
# implication to ETRS89 realization", EUREF2016<https://www.EUREF.EU/symposia/2016SanSebastian/01-02-Altamimi.pdf>}.
# _trfX(_ITRF2014_, _ETRF2000_, epoch=_F(2000),
#                             xform=_H(  53.7,     51.2,   -55.1,    1.02,      0.891,    5.39,     -8.712),
#                             rates=_H(   0.1,      0.1,    -1.9,    0.11,      0.081,    0.49,     -0.792))
_trfX(_ITRF2014_, _ETRF2000_, epoch=_F(2010),
                              xform=_H(  54.7,     52.2,   -74.1,    2.12,      1.701,   10.29,    -16.632),
                              rates=_H(   0.1,      0.1,    -1.9,    0.11,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2008_, _ETRF2000_, epoch=_F(2000),
#                             xform=_H(  52.1,     49.3,   -58.5,    1.34,      0.891,    5.39,     -8.712),
#                             rates=_H(   0.1,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2008_, _ETRF2000_, epoch=_F(2010),
                              xform=_H(  53.1,     50.3,   -76.5,    2.14,      1.701,   10.29,    -16.632),
                              rates=_H(   0.1,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2005_, _ETRF2000_, epoch=_F(2000),
#                             xform=_H(  54.1,     50.2,   -53.8,    0.4,       0.891,    5.39,     -8.712),
#                             rates=_H(  -0.2,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2005_, _ETRF2000_, epoch=_F(2010),
                              xform=_H(  52.1,     51.2,   -71.8,    1.2,       1.701,   10.29,    -16.632),
                              rates=_H(  -0.2,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2000_, _ETRF2000_, epoch=_F(2000),
#                             xform=_H(  54.0,     51.0,   -48.0,   _0_0,       0.891,    5.39,     -8.712),
#                             rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.792))
_trfX(_ITRF2000_, _ETRF2000_, epoch=_F(2010),
                              xform=_H(  54.0,     51.0,   -48.0,   _0_0,       1.701,   10.29,    -16.632),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.812))
_trfX(_ITRF97_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  54.0,     51.0,   -48.0,   -1.68,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF96_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  54.0,     51.0,   -48.0,   -1.68,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF94_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  47.3,     52.7,   -11.3,   -1.68,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF93_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H( 105.1,     48.9,   -13.9,   -2.17,      4.511,   13.67,    -17.032),
                              rates=_H(   2.9,      0.2,     0.6,   -0.01,      0.191,    0.68,     -0.862))
_trfX(_ITRF92_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  39.3,     50.7,    -3.3,   -0.97,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF91_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  27.3,     36.7,     2.7,   -2.37,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF90_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  29.3,     40.7,    18.7,   -2.67,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF89_,   _ETRF2000_, epoch=_F(2010),
                              xform=_H(  24.3,     16.7,    56.7,   -6.07,      1.701,   10.29,    -16.892),
                              rates=_H(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))

# GDA2020 "Geocentric Datum of Australia 2020 Technical Manual", v1.5, 2020-12-09, Table 3.3 and 3.4
# <https://www.ICSM.gov.AU/sites/default/files/2020-12/GDA2020%20Technical%20Manual%20V1.5_4.pdf>
# (the GDA2020 xforms are different but the rates are the same as GDA94, further below)
_trfX(_ITRF2014_, _GDA2020_,  epoch=_F(2020),
                              xform=_H(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,       1.50379,  1.18346,   1.20716))
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
                              xform=_H( 973.0,  -1919.2,  -482.9,   -0.9,     -25.79,    -9.65,    -11.66),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,      -0.053,    0.742,     0.032))
_trfX(_ITRF90_,   _WGS84_,    epoch=_F(1984),
                              xform=_H(  60.0,   -517.0,  -223.0,  -11.0,      18.3,     -0.3,       7.0),
                              rates=_H(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0))
del _H, _Hs

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
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
# date2epoch(2024,  1,  1) = 2024.003, epoch2date(2024.003) = 2024,  1,  1
# date2epoch(2024,  1, 15) = 2024.041, epoch2date(2024.041) = 2024,  1, 15
# date2epoch(2024,  1, 30) = 2024.082, epoch2date(2024.082) = 2024,  1, 30
# date2epoch(2024,  1, 31) = 2024.085, epoch2date(2024.085) = 2024,  1, 31
# date2epoch(2024,  2,  1) = 2024.087, epoch2date(2024.087) = 2024,  2,  2 *
# date2epoch(2024,  2, 15) = 2024.126, epoch2date(2024.126) = 2024,  2, 16 *
# date2epoch(2024,  2, 28) = 2024.161, epoch2date(2024.161) = 2024,  2, 28
# date2epoch(2024,  2, 29) = 2024.164, epoch2date(2024.164) = 2024,  3,  1 *
# date2epoch(2024,  3,  1) = 2024.167, epoch2date(2024.167) = 2024,  3,  2 *
# date2epoch(2024,  3, 15) = 2024.205, epoch2date(2024.205) = 2024,  3, 16 *
# date2epoch(2024,  3, 30) = 2024.246, epoch2date(2024.246) = 2024,  3, 31 *
# date2epoch(2024,  3, 31) = 2024.249, epoch2date(2024.249) = 2024,  4,  1 *
# date2epoch(2024,  4,  1) = 2024.251, epoch2date(2024.251) = 2024,  4,  1
# date2epoch(2024,  4, 15) = 2024.290, epoch2date(2024.290) = 2024,  4, 15
# date2epoch(2024,  4, 29) = 2024.328, epoch2date(2024.328) = 2024,  4, 29
# date2epoch(2024,  4, 30) = 2024.331, epoch2date(2024.331) = 2024,  4, 30
# date2epoch(2024,  5,  1) = 2024.333, epoch2date(2024.333) = 2024,  5,  1
# date2epoch(2024,  5, 15) = 2024.372, epoch2date(2024.372) = 2024,  5, 15
# date2epoch(2024,  5, 30) = 2024.413, epoch2date(2024.413) = 2024,  5, 30
# date2epoch(2024,  5, 31) = 2024.415, epoch2date(2024.415) = 2024,  6,  1 *
# date2epoch(2024,  6,  1) = 2024.418, epoch2date(2024.418) = 2024,  6,  2 *
# date2epoch(2024,  6, 15) = 2024.456, epoch2date(2024.456) = 2024,  6, 16 *
# date2epoch(2024,  6, 29) = 2024.495, epoch2date(2024.495) = 2024,  6, 30 *
# date2epoch(2024,  6, 30) = 2024.497, epoch2date(2024.497) = 2024,  7,  1 *
# date2epoch(2024,  7,  1) = 2024.500, epoch2date(2024.500) = 2024,  7,  1
# date2epoch(2024,  7, 15) = 2024.538, epoch2date(2024.538) = 2024,  7, 16 *
# date2epoch(2024,  7, 30) = 2024.579, epoch2date(2024.579) = 2024,  7, 30
# date2epoch(2024,  7, 31) = 2024.582, epoch2date(2024.582) = 2024,  7, 31
# date2epoch(2024,  8,  1) = 2024.585, epoch2date(2024.585) = 2024,  8,  1
# date2epoch(2024,  8, 15) = 2024.623, epoch2date(2024.623) = 2024,  8, 15
# date2epoch(2024,  8, 30) = 2024.664, epoch2date(2024.664) = 2024,  8, 31 *
# date2epoch(2024,  8, 31) = 2024.667, epoch2date(2024.667) = 2024,  9,  1 *
# date2epoch(2024,  9,  1) = 2024.669, epoch2date(2024.669) = 2024,  9,  2 *
# date2epoch(2024,  9, 15) = 2024.708, epoch2date(2024.708) = 2024,  9, 16 *
# date2epoch(2024,  9, 29) = 2024.746, epoch2date(2024.746) = 2024,  9, 30 *
# date2epoch(2024,  9, 30) = 2024.749, epoch2date(2024.749) = 2024, 10,  1 *
# date2epoch(2024, 10,  1) = 2024.751, epoch2date(2024.751) = 2024, 10,  1
# date2epoch(2024, 10, 15) = 2024.790, epoch2date(2024.790) = 2024, 10, 15
# date2epoch(2024, 10, 30) = 2024.831, epoch2date(2024.831) = 2024, 10, 30
# date2epoch(2024, 10, 31) = 2024.833, epoch2date(2024.833) = 2024, 10, 31
# date2epoch(2024, 11,  1) = 2024.836, epoch2date(2024.836) = 2024, 11,  1
# date2epoch(2024, 11, 15) = 2024.874, epoch2date(2024.874) = 2024, 11, 15
# date2epoch(2024, 11, 29) = 2024.913, epoch2date(2024.913) = 2024, 11, 29
# date2epoch(2024, 11, 30) = 2024.915, epoch2date(2024.915) = 2024, 12,  1 *
# date2epoch(2024, 12,  1) = 2024.918, epoch2date(2024.918) = 2024, 12,  2 *
# date2epoch(2024, 12, 15) = 2024.956, epoch2date(2024.956) = 2024, 12, 16 *
# date2epoch(2024, 12, 30) = 2024.997, epoch2date(2024.997) = 2024, 12, 31 *
# date2epoch(2024, 12, 31) = 2025.000, epoch2date(2025.000) = 2025,  1,  1 *
