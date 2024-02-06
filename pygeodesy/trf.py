
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

@var RefFrames.ETRF2000: RefFrame(name='ETRF2000', epoch=2005, datum=Datums.GRS80) .Xforms=(0, -14)
@var RefFrames.ETRF2005: RefFrame(name='ETRF2005', epoch=2005, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF2008: RefFrame(name='ETRF2008', epoch=2008, datum=Datums.GRS80) .Xforms=(0, 0)
@var RefFrames.ETRF2014: RefFrame(name='ETRF2014', epoch=2014, datum=Datums.GRS80) .Xforms=(0, -14)
@var RefFrames.ETRF2020: RefFrame(name='ETRF2020', epoch=2020, datum=Datums.GRS80) .Xforms=(0, -14)
@var RefFrames.ETRF88: RefFrame(name='ETRF88', epoch=1988, datum=Datums.GRS80) .Xforms=(0, 0)
@var RefFrames.ETRF89: RefFrame(name='ETRF89', epoch=1989, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF90: RefFrame(name='ETRF90', epoch=1990, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF91: RefFrame(name='ETRF91', epoch=1991, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF92: RefFrame(name='ETRF92', epoch=1992, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF93: RefFrame(name='ETRF93', epoch=1993, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF94: RefFrame(name='ETRF94', epoch=1994, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF96: RefFrame(name='ETRF96', epoch=1996, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.ETRF97: RefFrame(name='ETRF97', epoch=1997, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.GDA2020: RefFrame(name='GDA2020', epoch=2020, datum=Datums.GRS80) .Xforms=(0, -4)
@var RefFrames.GDA94: RefFrame(name='GDA94', epoch=1994, datum=Datums.GRS80) .Xforms=(0, -3)
@var RefFrames.ITRF2000: RefFrame(name='ITRF2000', epoch=1997, datum=Datums.GRS80) .Xforms=(15, -4)
@var RefFrames.ITRF2005: RefFrame(name='ITRF2005', epoch=2000, datum=Datums.GRS80) .Xforms=(8, -3)
@var RefFrames.ITRF2008: RefFrame(name='ITRF2008', epoch=2005, datum=Datums.GRS80) .Xforms=(17, -2)
@var RefFrames.ITRF2014: RefFrame(name='ITRF2014', epoch=2010, datum=Datums.GRS80) .Xforms=(16, -1)
@var RefFrames.ITRF2020: RefFrame(name='ITRF2020', epoch=2015, datum=Datums.GRS80) .Xforms=(16, 0)
@var RefFrames.ITRF88: RefFrame(name='ITRF88', epoch=1988, datum=Datums.GRS80) .Xforms=(3, -4)
@var RefFrames.ITRF89: RefFrame(name='ITRF89', epoch=1989, datum=Datums.GRS80) .Xforms=(4, -4)
@var RefFrames.ITRF90: RefFrame(name='ITRF90', epoch=1988, datum=Datums.GRS80) .Xforms=(6, -4)
@var RefFrames.ITRF91: RefFrame(name='ITRF91', epoch=1988, datum=Datums.GRS80) .Xforms=(4, -4)
@var RefFrames.ITRF92: RefFrame(name='ITRF92', epoch=1988, datum=Datums.GRS80) .Xforms=(4, -4)
@var RefFrames.ITRF93: RefFrame(name='ITRF93', epoch=1988, datum=Datums.GRS80) .Xforms=(4, -4)
@var RefFrames.ITRF94: RefFrame(name='ITRF94', epoch=1993, datum=Datums.GRS80) .Xforms=(4, -4)
@var RefFrames.ITRF96: RefFrame(name='ITRF96', epoch=1997, datum=Datums.GRS80) .Xforms=(5, -5)
@var RefFrames.ITRF97: RefFrame(name='ITRF97', epoch=1997, datum=Datums.GRS80) .Xforms=(5, -4)
@var RefFrames.NAD83: RefFrame(name='NAD83', epoch=1997, datum=Datums.GRS80) .Xforms=(0, -5)
@var RefFrames.WGS84: RefFrame(name='WGS84', epoch=1984, datum=Datums.GRS80) .Xforms=(0, -1)
@var RefFrames.WGS84g1150: RefFrame(name='WGS84g1150', epoch=2001, datum=Datums.GRS80) .Xforms=(0, 0)
@var RefFrames.WGS84g1674: RefFrame(name='WGS84g1674', epoch=2005, datum=Datums.GRS80) .Xforms=(0, 0)
@var RefFrames.WGS84g1762: RefFrame(name='WGS84g1762', epoch=2005, datum=Datums.GRS80) .Xforms=(0, 0)
'''

from pygeodesy.basics import map1, neg, isscalar, isstr, _xinstanceof, _zip
from pygeodesy.constants import _0_0s, _0_0, _0_001, _0_5, _1_0, _float as _F
from pygeodesy.datums import Datums, _earth_datum, _minus, Transform, _WGS84
from pygeodesy.errors import _IsnotError, TRFError, _xkwds
from pygeodesy.interns import MISSING, NN, _AT_, _COMMASPACE_, _cartesian_, _conversion_, \
                             _datum_, _DOT_, _ellipsoidal_, _exists_, _invalid_, _MINUS_, \
                             _NAD83_, _no_, _PLUS_, _reframe_, _s_, _scalar_, _SPACE_, \
                             _sx_, _sy_, _sz_, _to_, _tx_, _ty_, _tz_, _WGS84_, _x_, \
                             _intern as _i
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import ADict, classname, _lazyNamedEnumItem as _lazy, _Named, \
                           _NamedEnum, _NamedEnumItem, _NamedTuple,  Fmt, unstr
from pygeodesy.props import Property, Property_RO, property_RO
# from pygeodesy.streprs import Fmt, unstr  # from .named
from pygeodesy.units import Epoch, Float

from math import ceil

__all__ = _ALL_LAZY.trf
__version__ = '24.02.04'

_EP0CH =  Epoch(0, low=0)
_Es    = {_EP0CH: _EP0CH}  # L{Epoch} de_dup
_Ps    = {}  # L{TRFXform7Tuple} de-dup, deleted below

_GRS80    =  Datums.GRS80
_inverse_ = 'inverse'
_MAS      = _MM = _PPB = Float  # Units
_MM2M     = _0_001  # scale mm2m, ppb2ppM, mas2as
_mDays    = (0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0)
_366_0    = _F(sum(_mDays))
_v_       = 'v'

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
_ETRF2020_   = _i('ETRF2020')
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


def _E(epoch):
    '''(INTERNAL) De-dup L{Epochs}s.
    '''
    e = Epoch(_F(epoch))
    return _Es.setdefault(e, e)


class RefFrame(_NamedEnumItem):
    '''Terrestrial Reference Frame (TRF) parameters.
    '''
    _datum = _GRS80  # Datums.GRS80 or .WGS84 (L{Datum})
    _epoch = _EP0CH  # epoch, fractional year (L{Epoch})

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
        self._Xs = {}
        if name:
            self._register(RefFrames, name)

    def __eq__(self, other):
        return isinstance(other, RefFrame) and other.epoch == self.epoch \
                                           and other.datum == self.datum

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

    def toRefFrame(self, point, reframe2, **epoch2_epoch_name):
        '''Convert a cartesian or geodetic point from this to another reframe and epoch.

           @return: A copy of the B{C{point}}, converted or renamed.

           @see: Ellipsoidal methods L{LatLon.toRefFrame<ellipsoidalBase.LatLonEllipsoidalBase.toRefFrame>}
                 and L{Cartesian.toRefFrame<ellipsoidalBase.CartesianEllipsoidalBase.toRefFrame>}
                 for more details.
        '''
        b = _MODS.ellipsoidalBase
        _xinstanceof(b.LatLonEllipsoidalBase, b.CartesianEllipsoidalBase, point=point)
        p = point.dup(reframe=self) if point.reframe != self else point
        return p.toRefFrame(reframe2, **_xkwds(epoch2_epoch_name, name=self.name))

    def toStr(self, epoch=None, name=NN, **unused):  # PYCHOK expected
        '''Return this reference frame as a text string.

           @kwarg epoch: Override this reframe's epoch (C{scalar} or C{str}).
           @kwarg name: Override name (C{str}) or C{None} to exclude the
                        reframe's name.

           @return: This L{RefFrame}'s attributes (C{str}).
        '''
        D = self.datum
        e = self.epoch if epoch is None else _Epoch(epoch)
        t = (Fmt.EQUAL(name=repr(name or self.name)),
             Fmt.EQUAL(epoch=e),
             Fmt.EQUAL(datum=_DOT_(classname(D) + _s_, D.name)))
        return _COMMASPACE_.join(t[1:] if name is None else t)

    def Xform(self, reframe2):
        '''Get the converter Xform I{from} this reference frame I{to} C{reframe2}.

           @arg reframe2: Destination frame to convert I{to} (L{RefFrame} or C{str}).

           @return: The L{TRFXform} instance or C{None} if not available.

           @raise TypeError: Invalid B{C{reframe2}}.
        '''
        _xinstanceof(RefFrame, str, reframe2=reframe2)
        n2 = reframe2 if isstr(reframe2) else reframe2.name
        return self._Xs.get(n2, None)

    def Xforms(self, inverse=False):
        '''Return all Xforms converting I{from} or I{to} this reference frame.

           @kwarg inverse: If C{True}, get all I{inverse, to} Xforms (C{bool}).

           @return: An L{ADict} of I{from} C{[refName2]=}L{TRFXform}s) or if
                    C{B{inverse}=True} of I{to} C{[refName1]=}L{TRFXform}s).
        '''
        def _Xi(n2):
            for n1, r in RefFrames.items():
                if n2 in r._Xs:
                    yield n1, r

        return ADict(_Xi(self.name)) if inverse else ADict(self._Xs)


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
    ETRF88     = _lazy(_ETRF88_,     _E(1988)),  # epoch, datum?
    ETRF89     = _lazy(_ETRF89_,     _E(1989)),  # epoch, datum?
    ETRF90     = _lazy(_ETRF90_,     _E(1990)),  # epoch, datum?
    ETRF91     = _lazy(_ETRF91_,     _E(1991)),  # epoch, datum?
    ETRF92     = _lazy(_ETRF92_,     _E(1992)),  # epoch, datum?
    ETRF93     = _lazy(_ETRF93_,     _E(1993)),  # epoch, datum?
    ETRF94     = _lazy(_ETRF94_,     _E(1994)),  # epoch, datum?
    ETRF96     = _lazy(_ETRF96_,     _E(1996)),  # epoch, datum?
    ETRF97     = _lazy(_ETRF97_,     _E(1997)),  # epoch, datum?
    ETRF2000   = _lazy(_ETRF2000_,   _E(2005)),
    ETRF2005   = _lazy(_ETRF2005_,   _E(2005)),  # epoch, datum?
    ETRF2008   = _lazy(_ETRF2008_,   _E(2008)),  # epoch, datum?
    ETRF2014   = _lazy(_ETRF2014_,   _E(2014)),  # epoch, datum?
    ETRF2020   = _lazy(_ETRF2020_,   _E(2020)),  # epoch, datum?
    GDA94      = _lazy(_GDA94_,      _E(1994)),  # Australia
    GDA2020    = _lazy(_GDA2020_,    _E(2020)),  # Australia
    ITRF88     = _lazy(_ITRF88_,     _E(1988)),
    ITRF89     = _lazy(_ITRF89_,     _E(1989)),
    ITRF90     = _lazy(_ITRF90_,     _E(1988)),
    ITRF91     = _lazy(_ITRF91_,     _E(1988)),
    ITRF92     = _lazy(_ITRF92_,     _E(1988)),
    ITRF93     = _lazy(_ITRF93_,     _E(1988)),
    ITRF94     = _lazy(_ITRF94_,     _E(1993)),
    ITRF96     = _lazy(_ITRF96_,     _E(1997)),
    ITRF97     = _lazy(_ITRF97_,     _E(1997)),
    ITRF2000   = _lazy(_ITRF2000_,   _E(1997)),  # aka ITRF00
    ITRF2005   = _lazy(_ITRF2005_,   _E(2000)),
    ITRF2008   = _lazy(_ITRF2008_,   _E(2005)),  # aka ITRF08
    ITRF2014   = _lazy(_ITRF2014_,   _E(2010)),
    ITRF2020   = _lazy(_ITRF2020_,   _E(2015)),
    NAD83      = _lazy(_NAD83_,      _E(1997)),  # aka CORS96
    WGS84      = _lazy(_WGS84_,      _E(1984), _WGS84),
    WGS84g1150 = _lazy(_WGS84g1150_, _E(2001), _WGS84),
    WGS84g1674 = _lazy(_WGS84g1674_, _E(2005), _WGS84),
    WGS84g1762 = _lazy(_WGS84g1762_, _E(2005), _WGS84))  # same epoch


class TransformXform(Transform):
    '''Helmert transformation, extended with an C{Xform} converter.

       @see: L{Transform<datums.Transform>} and L{Xform<TRFXform>}.
    '''
    _Xform = None

    def __init__(self, name=NN, **tx_ty_tz_s_sx_sy_sz):  # PYCHOK signature
        '''New L{TransformXform}.

           @kwarg name: Optional name (C{str}), I{not registered}.

           @see: L{Transform<datums.Transform>} for details.

           @note: The name of this L{TransformXform} starts with C{"-"}
                  if its C{Xform} was applied I{inversed}.
        '''
        Transform.__init__(self, **tx_ty_tz_s_sx_sy_sz)
        if name:
            self.name = name

#   def _rps2(self, s_):
#       '''(INTERNAL) Rate to rate.
#       '''
#       return s_, s_

    def toRefFrame(self, point, epoch=None, **epoch2_datum_name):  # PYCHOK signature
        '''Convert a cartesian or geodetic point using this transform and C{refName1},
           C{refName2} and C{epoch} of this transform's C{Xform}.

           @kwarg epoch: Optional epoch ((L{Epoch}, C{scalar} or C{str})), overriding
                         this transform's Xform's C{epoch}.

           @return: A copy of the B{C{point}}, converted or renamed.

           @see: Method L{RefFrame.toRefFrame<trf.RefFrame.toRefFrame>} for more details.
        '''
        return self.Xform.toRefFrame(point, epoch=epoch, **_xkwds(epoch2_datum_name))

    def velocities(self, factor=_MM2M):
        '''Compute the X, Y and Z I{velocities} of this transform.

           @kwarg factor: Factor to scale this Xform's C{rates} (C{scalar}), default
                          from C{milli-meter-} to C{meter-per-year}.

           @return: A L{Vector3Tuple}C{(x, y, z)} or C{None}.

           @raise TypeError: Invalid B{C{factor}}.

           @see: Alamimi, Z. "EUREF-TN-1-Jan-31-2024", U{Appendix A, equation (3)
                 <http://ETRS89.ENSG.IGN.FR/pub/EUREF-TN-1-Jan-31-2024.pdf>}.
        '''
        v = self.Xform
        if v is not None:
            r = v.rates * factor  # eq (3) ...
            T = self.dup(tx=r.sx, ty=r.sy, tz=r.sz,  # Xyy-dot?
                         s1=0, name=NN(self.name, _v_))
            v = T.transform(r.tx, r.ty, r.tz)  # Xyy?
        return v

    @Property
    def Xform(self):
        '''Get the Xform of this Helmert (L{TRFXform} or C{None}).
        '''
        return self._Xform

    @Xform.setter  # PYCHOK setter!
    def Xform(self, Xform):
        '''Set the Xform of this Helmert (L{TRFXform}).

           @raise TypeError: Invalid B{C{Xform}}.
        '''
        _xinstanceof(TRFXform, Xform=Xform)
        self._Xform = Xform


class TRFXform7Tuple(_NamedTuple):
    '''7-Tuple C{(tx, ty, tz, s, sx, sy, sz)} of conversion parameters with
       translations C{tx}, C{ty} and C{tz} in C{milli-meter}, scale C{s} in
       C{ppb} and rates C{sx}, C{sy} and C{sz} in C{milli-arc-seconds-per-year}.

       @note: The parameters are also named as C{(Tx, Ty, Tz, D, Rx, Ry, Rz)}
              or C{(T1, T2, T3, D, R1, R2, R3)}.

       @see: Class L{TransformXform}'s matching keyword argument names.
    '''
    _Names_ = (_tx_, _ty_, _tz_, _s_,  _sx_, _sy_, _sz_)  # == kwds
    _Units_ = (_MM,  _MM,  _MM,  _PPB, _MAS, _MAS, _MAS)

    def __add__(self, other):
        '''Return the sum of this and an other L{TRFXform7Tuple}.
        '''
        _xinstanceof(TRFXform7Tuple, other=other)
        return TRFXform7Tuple(((s + t) for s, t in _zip(self, other)),
                              name=_PLUS_)  # .fsums._add_op_

    def __mul__(self, factor):
        if not isscalar(factor):
            raise _IsnotError(_scalar_, factor=factor)
        return type(self)(_ * factor for _ in self)

    def __neg__(self):
        return self.inverse()

    def inverse(self, name=NN):
        '''Return the inverse of this transform.

           @kwarg name: Optional name (C{str}).

           @return: Inverse (L{TransformXform}).
        '''
        n = name or _minus(self.name)
        return type(self)(map(neg, self), name=n)

    @Property_RO
    def isunity(self):
        '''Is this a C{unity, identidy} transform (C{bool}).
        '''
        return not any(self)


class TRFXform(_Named):
    '''Terrestrial Reference Frame (TRF) converter between two
       reference frames observed at an C{epoch}.
    '''
    epoch    = _EP0CH
    refName1 =  \
    refName2 =  NN
    rates    =  \
    xform    =  TRFXform7Tuple(_0_0s(len(
                TRFXform7Tuple._Names_)))

    def __init__(self, refName1, refName2, epoch=None,
                                           xform=None, rates=None, name=NN):
        '''New L{trfXform}.

           @arg refName1: Source reframe (C{str}), converting I{from}.
           @arg refName2: Destination reframe (C{str}), converting I{to}.
           @kwarg epoch: Epoch, a fractional year (L{Epoch}, C{scalar} or C{str}).
           @kwarg xform: I{Transform} parameters (L{TRFXform7Tuple}).
           @kwarg rates: I{Rate} parameters (L{TRFXform7Tuple}), as B{C{xform}},
                         but in C{units-per-year}.
           @kwarg name: Optional name (C{str})

           @return: The new TRF converter (L{TRFXform}).

           @raise TRFError: Invalid B{C{epoch}}, B{C{xform}} or B{C{rates}}.
        '''
        self.refName1 = str(refName1)
        self.refName2 = str(refName2)
        _xinstanceof(TRFXform7Tuple, xform=xform, rates=rates)
        self.epoch = _Epoch(epoch)
        self.xform =  xform
        self.rates =  rates
        self.name  =  name or self.toStr()

    def __lt__(self, other):
        _xinstanceof(TRFXform, other=other)
        return self.refName1 < other.refName1 or \
               self.refName2 < other.refName2 or \
               self.epoch    < other.epoch

    def __neg__(self):
        return self.inverse()

    def __repr__(self):
        return self.toRepr()

    def _at(self, epoch):
        '''(INTERNAL) Return C{self.name@epoch}.
        '''
        n = self.name
        if epoch != self.epoch:
            _t = _AT_(NN, epoch)
            if not n.endswith(_t):
                n = NN(n, _t)
        return n

    def _7Helmert(self, epoch, mm2m):
        '''(INTERNAL) Get the named Helmert parameters observed
           at C{epoch} converted by C{mm2m}, I{signed}.
        '''
        e  = epoch - self.epoch  # delta in fractional years
        ns = TRFXform7Tuple._Names_
        return dict((n, (x + e * r) * mm2m) for n, x, r in
                    _zip(ns, self.xform, self.rates))  # strict=True

    def inverse(self, name=NN):
        '''Return the inverse of this Xform.

           @kwarg name: Optional name (C{str}).

           @return: Inverse (L{TRFXform}).
        '''
        return self.dup(refName1=self.refName2, xform=-self.xform,
                        refName2=self.refName1, rates=-self.rates,
                        name=name or _minus(self.name))

    def _range(self, *e1_e2):
        '''(INTERNAL) Epoch "range", distance.
        '''
        def _r(e):
            e -= self.epoch
            return (-e) if e < 0 else (e * 2)

        return sum(map(_r, e1_e2)) if e1_e2 else self.epoch

    def toEpoch(self, epoch):
        '''Convert this Xform to B{C{epoch}}, if needed.

           @arg epoch: Epoch, fractional year (C{Epoch}, C{scalar}
                       or C{str}).

           @return: This Xform or a copy converted to B{C{epoch}}.

           @raise TRFError: Invalid B{C{epoch}}.
        '''
        e = _Epoch(epoch) - self.epoch
        if e:
            t = TRFXform7Tuple((x + e * r) for x, r in
                               _zip(self.xform, self.rates))
            X = self.dup(epoch=epoch, xform=t, name=self._at(epoch))
        else:
            X = self
        return X

    def toRefFrame(self, point, epoch=None, datum=_GRS80, **epoch2_name):
        '''Convert a cartesian or geodetic point using this Xform's C{refName1},
           C{refName2}, C{epoch} and L{Transform<trf.TransformXform>}.

           @kwarg epoch: Optional epoch ((L{Epoch}, C{scalar} or C{str})), overriding
           this transform's Xform's C{epoch}.

           @return: A copy of the B{C{point}}, converted or renamed.

           @see: Method L{RefFrame.toRefFrame<trf.RefFrame.toRefFrame>} for more details.
        '''
        def _r(name):
            r = RefFrames.get(name)
            if r is None or r.datum != datum:
                r = RefFrame(self.epoch, datum)
                r.name = name  # unregistered
            return r

        return _r(self.refName1).toRefFrame(point,
               _r(self.refName2), epoch=epoch, **_xkwds(epoch2_name))

    def toRepr(self, **unused):  # PYCHOK signature
        '''Return the represention of this Xform (C{str}).
        '''
        return unstr(self.classname, epoch=self.epoch, name=self.name)

    def toStr(self, epoch=None, **unused):  # PYCHOK signature
        '''Return this Xform as C{"fromRef@epoch_toRef"} (C{str}).

           @kwarg epoch: Epoch (C{scalar}), overriding this Xform's epoch.
        '''
        e = epoch or self.epoch
        return NN(_AT_(self.refName1, e), _x_, self.refName2)

    def toTransform(self, epoch1, epoch2=None, inverse=False):
        '''Combine this Xform observed at B{C{epoch1}} into a Helmert
           L{TransformXform}, optionally at B{C{epoch2 or epoch1}}.

           @arg epoch1: Epoch to observe I{from} (C{scalar}).
           @kwarg epoch2: Optional epoch to observe I{to} (C{scalar}).
           @kwarg inverse: Invert the Helmert tranformation (C{bool}).

           @raise TRFError: Invalid B{C{epoch1}} or B{C{epoch2}}.

           @note: Translations are converted from C{milli-meter} to
                  C{meter}, shift rates from C{milli-arc-seconds} to
                  C{arc-seconds} and scale from C{ppb} to C{ppM}.
        '''
        e = Epoch(epoch2=epoch1 if epoch2 is None else epoch2)
        X = self if e == epoch1 else self.toEpoch(epoch1)
        if inverse:
            X = X.inverse()
        H = TransformXform(name=X._at(e), **X._7Helmert(e, _MM2M))
        H.Xform = X
        return H


def date2epoch(year, month, day):
    '''Return the C{epoch} for a calendar day.

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
            e = y + float(sum(_mDays[:m]) + d) / _366_0
            return Epoch(e, low=0)

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
    e = _Epoch(epoch, low=0)
    y =  int(e)
    d =  int(ceil(_366_0 * (e - y)))
    for m, n in enumerate(_mDays[1:]):
        if d > n:
            d -= n
        else:
            break
    return y, (m + 1), max(1, d)


def _Epoch(e, **kwds):
    '''(INTERNAL) Get an C{Epoch(B{e})}.
    '''
    return Epoch(e, **kwds) if kwds or not isinstance(e, Epoch) else e


def _eT0Ds4(inst, reframe, epoch, reframe2, epoch2):
    '''(INTERNAL) Get epoch, a Helmert L{Transform}s or 0-tuple, datum
       and datum2 to convert B{C{refFrame}} observed at B{C{epoch}}
       into B{C{refFrame2}} observed at C{epoch2 or epoch}.
    '''
    r = reframe or inst.reframe
    if not r:
        t = _SPACE_(_DOT_(repr(inst), _reframe_), MISSING)
        raise TRFError(_no_(_conversion_), txt=t)

    _xinstanceof(RefFrame, reframe2=reframe2, reframe=r)

    e1 = _Epoch(epoch or inst.epoch or r.epoch)
    e2 =  e1 if epoch2 is None else Epoch(epoch2=epoch2)

    t0 = _toTransform0(r.name, e1, reframe2.name, e2)
    if t0 is None:
        t = _SPACE_(RefFrame.__name__, _AT_(r.name, e1),
                          _to_, _AT_(reframe2.name, e2))
        raise TRFError(_no_(_conversion_), txt=t)

    return e2, t0, r.datum, reframe2.datum


def _immediate(n1, n2, *unused):
    '''(INTERNAL) Get the C{n1} to C{n2} Xform or C{None}.
    '''
    X = RefFrames.get(n1)
    if X:
        X = X._Xs.get(n2, None)
    return X


def _intermediate(n1, n2, *e1_e2):
    '''(INTERNAL) Find an Xform "between" C{n1} and C{n2}.
    '''
    def _e(e1=0, e2=0):
        return e1 - e2

    def _X4(_Xs):
        for n, X1 in _Xs.items():
            X2 = _immediate(n, n2)
            if X2:
                yield (X1._range(*e1_e2) + X2._range(*e1_e2)), X1, X2, True
            X2 = _immediate(n2, n)
            if X2:
                yield (X1._range(*e1_e2) + X2._range(*e1_e2)), X1, X2, False

    X, r = None, RefFrames.get(n1)
    if r:
        i = sorted(_X4(r._Xs))
        if i:
            _, X1, X2, f = i.pop(0)
            e1, e2 = X1.epoch, X2.epoch
            if e1 > e2 and _e(*e1_e2) < 0:  # X1 to e2
                X1 = X1.toEpoch(e2)
                e1 = e2
            elif e1 != e2:  # X2 to e1
                X2 = X2.toEpoch(e1)
            # <https://Geodesy.NOAA.gov/TOOLS/Htdp/Pearson_Snay_2012.pdf,
            # Apendix, Table 7, last column "Sum of the previous ..."
            X = TRFXform(n1, n2, epoch=e1,
                                 xform=X1.xform + (X2.xform if f else -X2.xform),
                                 rates=X1.rates + (X2.rates if f else -X2.rates),
                                 name=NN(X1._at(e1), _PLUS_ if f else _MINUS_,
                                         X2._at(e1)))
    return X


def _P(*ps):  # deleted below
    '''(INTERNAL) De-dup L{TRFXform7Tuple}s.
    '''
    t = TRFXform7Tuple(map(_F, ps))
    return _Ps.setdefault(t, t)  # PYCHOK del


def _reframe(**name_reframe):
    '''(INTERNAL) Get a C{reframe}.
    '''
    if len(name_reframe) == 1:
        for n, r in name_reframe.items():
            break
        if isstr(r):
            r = RefFrames.get(r)
        if r and isinstance(r, RefFrame):
            return r
    raise TRFError(**name_reframe)  # _invalid_


def _toTransform0(n1, e1, n2, e2):
    '''(INTERNAL) Return a L{TransformXform}, 0-tuple or C{None}.
    '''
    if n1 == n2 or (n1.startswith(_ITRF_) and n2.startswith(_WGS84_)) \
                or (n2.startswith(_ITRF_) and n1.startswith(_WGS84_)):
        return ()

    for _f in (_immediate, _intermediate):
        X = _f(n1, n2, e1, e2)
        Y = _f(n2, n1, e1, e2)
        if X and Y:
            if X._range(e1, e2) < Y._range(e1, e2):
                Y = None
            else:
                X = None
        if X:
            return X.toTransform(e1, epoch2=e2)
        if Y:
            return Y.toTransform(e1, epoch2=e2, inverse=True)

    return None


def trfTransform0(reframe, reframe2, epoch=None, epoch2=None):
    '''Get the L{TransformXform} to convert a C{reframe} observed at C{epoch} to
       C{reframe2} at observed at C{epoch2 or epoch}.

       @arg reframe: Reference frame to convert I{from} (L{RefFrame} or C{str}).
       @arg reframe2: Reference frame to convert I{to} (L{RefFrame} or C{str}).
       @arg epoch: Epoch to observe I{from} (L{Epoch}, C{scalar} or C{str}),
                   otherwise C{B{reframe}}'s C{epoch}.
       @kwarg epoch2: Optional epoch to observe to observe I{to} (L{Epoch}, C{scalar}
                      or C{str}), otherwise B{C{epoch}} or C{B{reframe}}'s C{epoch}.

       @return: A L{TransformXform} instance, a C{0-tuple} for I{unity, identity} or
                C{None} if no conversion exists.

       @raise TRFError: Invalid B{C{reframe}}, B{C{epoch}}, B{C{reframe2}} or
                        B{C{epoch2}}.

       @raise TypeError: Invalid B{C{reframe}} or B{C{reframe2}}.
    '''
    r1 = _reframe(reframe=reframe)
    e1 =  r1.epoch if epoch is None else _Epoch(epoch)
    r2 = _reframe(reframe2=reframe2)
    e2 =  e1 if epoch2 is None else Epoch(epoch2=epoch2)
    return _toTransform0(r1.name, e1, r2.name, e2)


def trfXform(reframe1, reframe2, epoch=None, xform=None, rates=None):
    '''Define a new Terrestrial Reference Frame (TRF) converter.

       @arg reframe1: Source reframe (L{RefFrame} or C{str}), converting I{from}.
       @arg reframe2: Destination reframe (L{RefFrame} or C{str}), converting I{to}.
       @kwarg epoch: Epoch, a fractional year (L{Epoch}, C{scalar} or C{str})
                     or C{None} for C{B{reframe2}}'s epoch.
       @kwarg xform: I{Transform} parameters (L{TRFXform7Tuple}).
       @kwarg rates: I{Rate} parameters (L{TRFXform7Tuple}), as B{C{xform}},
                     but in C{units-per-year}.

       @return: The new TRF converter (L{TRFXform}).

       @raise TRFError: Invalid B{C{reframe1}}, B{C{reframe2}}, B{C{epoch}},
                        B{C{xform}} or B{C{rates}} or the TRF converter
                        already exists.
    '''
    r1 = _reframe(reframe1=reframe1)
    r2 = _reframe(reframe2=reframe2)
    e  =  r2.epoch if epoch is None else _Epoch(epoch)
    _xinstanceof(TRFXform7Tuple, xform=xform, rates=rates)
    return _trfX(r1.name, r2.name, epoch=e, xform=xform, rates=rates)


def _trfX(n1, n2, **epoch_xform_rates):
    '''(INTERNAL) New, I{unique} L{TRFXform} converter.
    '''
    r1 = RefFrames.get(n1)
    r2 = RefFrames.get(n2)
    if r1 and r2:
        if n2 in r1._Xs:
            t = _exists_
        elif n1 in r2._Xs:
            t = _SPACE_(_inverse_, _exists_)
        else:
            r1._Xs[n2] = X = TRFXform(n1, n2, **epoch_xform_rates)
            return X
    else:
        t = _invalid_
    raise TRFError(trfX=_SPACE_(n1, _to_, n2), txt=t)  # _NameError


# TRF conversions specified as an epoch and 2 sets of 7-parameters.  Most from Altamimi, Z. U{"EUREF Technical
# Note 1: Relationship and Transformation between the International and the European Terrestrial Reference
# Systems"<https://ERTS89.ENSG.IFN.Fr/pub/EUREF-TN-1-Jan-31-2024.pdf>} Appendix A, more at U{Quinsy QPS <https://
# confluence.QPS.NL/qinsy/files/latest/en/182618383/182618384/1/1579182881000/ITRF_Transformation_Parameters.xlsx>}.
# See also U{Quinsy International Terrestrial Reference Frame 2014 (ITRF2014)<https://confluence.QPS.NL/qinsy/latest/
# en/international-terrestrial-reference-frame-2014-itrf2014-182618383.html>}.
_trfX(_ITRF2020_, _ITRF2014_, epoch=_E(2015),  # <https://ITRF.IGN.Fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt>
                              xform=_P(  -1.4,     -0.9,     1.4,   -0.42,     _0_0,      _0_0,     _0_0),
                              rates=_P(  _0_0,     -0.1,     0.2,    _0_0,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2008_, epoch=_E(2015),
                              xform=_P(   0.2,      1.0,     3.3,   -0.29,     _0_0,      _0_0,     _0_0),
                              rates=_P(  _0_0,     -0.1,     0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2005_, epoch=_E(2015),
                              xform=_P(   2.7,      0.1,    -1.4,    0.65,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.3,     -0.1,     0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF2000_, epoch=_E(2015),
                              xform=_P(  -0.2,      0.8,   -34.2,    2.25,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.1,     _0_0,    -1.7,    0.11,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2020_, _ITRF97_,   epoch=_E(2015),
                              xform=_P(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF96_,   epoch=_E(2015),
                              xform=_P(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF94_,   epoch=_E(2015),
                              xform=_P(   6.5,     -3.9,   -77.9,    3.98,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF93_,   epoch=_E(2015),
                              xform=_P( -65.8,      1.9,   -71.3,    4.47,     -3.36,     -4.33,     0.75),
                              rates=_P(  -2.8,     -0.2,    -2.3,    0.12,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2020_, _ITRF92_,   epoch=_E(2015),
                              xform=_P(  14.5,     -1.9,   -85.9,    3.27,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF91_,   epoch=_E(2015),
                              xform=_P(  26.5,     12.1,   -91.9,    4.67,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF90_,   epoch=_E(2015),
                              xform=_P(  24.5,      8.1,  -107.9,    4.97,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF89_,   epoch=_E(2015),
                              xform=_P(  29.5,     32.1,  -145.9,    8.37,     _0_0,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2020_, _ITRF88_,   epoch=_E(2015),
                              xform=_P(  24.5,     -3.9,  -169.9,   11.47,      0.1,      _0_0,      0.36),
                              rates=_P(   0.1,     -0.6,    -3.1,    0.12,     _0_0,      _0_0,      0.02))

# see U{Transformation Parameters ITRF2014<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>} and
# Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Appendix A.
_trfX(_ITRF2014_, _ITRF2008_, epoch=(2010),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2014/tp_14-08.php>
                              xform=_P(   1.6,      1.9,     2.4,   -0.02,     _0_0,      _0_0,     _0_0),
                              rates=_P(  _0_0,     _0_0,    -0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF2005_, epoch=_E(2010),
                              xform=_P(   2.6,     _1_0,    -2.3,    0.92,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.3,     _0_0,    -0.1,    0.03,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF2000_, epoch=_E(2010),
                              xform=_P(   0.7,      1.2,   -26.1,    2.12,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.1,      0.1,    -1.9,    0.11,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2014_, _ITRF97_,   epoch=_E(2010),
                              xform=_P(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF96_,   epoch=_E(2010),
                              xform=_P(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF94_,   epoch=_E(2010),
                              xform=_P(   7.4,     -0.5,   -62.8,    3.8,      _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF93_,   epoch=_E(2010),
                              xform=_P( -50.4,      3.3,   -60.2,    4.29,     -2.81,     -3.38,     0.4),
                              rates=_P(  -2.8,     -0.1,    -2.5,    0.12,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2014_, _ITRF92_,   epoch=_E(2010),
                              xform=_P(  15.4,      1.5,   -70.8,    3.09,     _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF91_,   epoch=_E(2010),
                              xform=_P(  27.4,     15.5,   -76.8,    4.49,     _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF90_,   epoch=_E(2010),
                              xform=_P(  25.4,     11.5,   -92.8,    4.79,     _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF89_,   epoch=_E(2010),
                              xform=_P(  30.4,     35.5,  -130.8,    8.19,     _0_0,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2014_, _ITRF88_,   epoch=_E(2010),
                              xform=_P(  25.4,     -0.5,  -154.8,   11.29,      0.1,      _0_0,      0.26),
                              rates=_P(   0.1,     -0.5,    -3.3,    0.12,     _0_0,      _0_0,      0.02))

# Pearson, C. & Snay, R. U{"Introducing HTDP 3.1 to transform coordinates across time and spatial reference
# frames"<https://Geodesy.NOAA.gov/TOOLS/Htdp/Pearson_Snay_2012.pdf> Table 7, 1st and 2nd column
_trfX(_ITRF2008_, _ITRF2005_, epoch=_E(2005),
                              xform=_P(  -0.5,     -0.9,    -4.7,    0.94,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,     _0_0))
# _trfX(_ITRF2008_, _ITRF2005_, epoch=_E(1997),
#                             xform=_P(  -2.9,     -0.9,    -4.7,    0.94,     _0_0,      _0_0,     _0_0),
#                             rates=_P(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,     _0_0))
# see U{Transformation Parameters ITRF2008<http://ITRF.IGN.Fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt>}
# _trfX(_ITRF2008_, _ITRF2005_, epoch=_E(2000),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2008/tp_08-05.php>
#                             xform=_P(  -2.0,     -0.9,    -4.7,    0.94,     _0_0,      _0_0,     _0_0),
#                             rates=_P(   0.3,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,     _0_0))
_trfX(_ITRF2008_, _ITRF2000_, epoch=_E(2000),
                              xform=_P(  -1.9,     -1.7,   -10.5,    1.34,     _0_0,      _0_0,     _0_0),
                              rates=_P(   0.1,      0.1,    -1.8,    0.08,     _0_0,      _0_0,     _0_0))
_trfX(_ITRF2008_, _ITRF97_,   epoch=_E(2000),
                              xform=_P(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF96_,   epoch=_E(2000),
                              xform=_P(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF94_,   epoch=_E(2000),
                              xform=_P(   4.8,      2.6,   -33.2,    2.92,     _0_0,      _0_0,      0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF93_,   epoch=_E(2000),
                              xform=_P( -24.0,      2.4,   -38.6,    3.41,     -1.71,     -1.48,    -0.3),
                              rates=_P(  -2.8,     -0.1,    -2.4,    0.09,     -0.11,     -0.19,     0.07))
_trfX(_ITRF2008_, _ITRF92_,   epoch=_E(2000),
                              xform=_P(  12.8,      4.6,   -41.2,    2.21,     _0_0,      _0_0,      0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,      _0_0,      0.02))
_trfX(_ITRF2008_, _ITRF91_,   epoch=_E(2000),
                              xform=_P(  24.8,     18.6,   -47.2,    3.61,     _0_0,     _0_0,       0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF90_,   epoch=_E(2000),
                              xform=_P(  22.8,     14.6,   -63.2,    3.91,     _0_0,     _0_0,       0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF89_,   epoch=_E(2000),
                              xform=_P(  27.8,     38.6,  -101.2,    7.31,     _0_0,     _0_0,       0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2008_, _ITRF88_,   epoch=_E(2000),
                              xform=_P(  22.8,      2.6,  -125.2,   10.41,      0.1,     _0_0,       0.06),
                              rates=_P(   0.1,     -0.5,    -3.2,    0.09,     _0_0,     _0_0,       0.02))

# Pearson, C. & Snay, R. U{"Introducing HTDP 3.1 to transform coordinates across time and spatial reference
# frames"<https://Geodesy.NOAA.gov/TOOLS/Htdp/Pearson_Snay_2012.pdf> Table 7, 3rd column
# _trfX(_ITRF2005_, _ITRF2000_, epoch=_E(1997),
#                             xform=_P(   0.7,     -1.1,    -0.4,    0.16,     -0.2,      0.1       -1.8),
#                             rates=_P(  -0.2,      0.1,    -1.8,    0.08,     _0_0,     _0_0,      _0_0))
_trfX(_ITRF2005_, _ITRF2000_, epoch=_E(2000),  # <http://ITRF.ENSG.IGN.Fr/ITRF_solutions/2005/tp_05-00.php>
                              xform=_P(   0.1,     -0.8,    -5.8,    0.4,      _0_0,     _0_0,      _0_0),
                              rates=_P(  -0.2,      0.1,    -1.8,    0.08,     _0_0,     _0_0,      _0_0))

_trfX(_ITRF2000_, _ITRF97_,   epoch=_E(1997),
                              xform=_P(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,      -0.02))  # 0.02?
_trfX(_ITRF2000_, _ITRF96_,   epoch=_E(1997),
                              xform=_P(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF94_,   epoch=_E(1997),
                              xform=_P(   0.67,     0.61,   -1.85,   1.55,     _0_0,     _0_0,      _0_0),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF93_,   epoch=_E(1988),
                              xform=_P(  12.7,      6.5,   -20.9,    1.95,     -0.39,     0.8,      -1.14),
                              rates=_P(  -2.9,     -0.2,    -0.6,    0.01,     -0.11,    -0.19,      0.07))
_trfX(_ITRF2000_, _ITRF92_,   epoch=_E(1988),
                              xform=_P(   1.47,     1.35,   -1.39,   0.75,     _0_0,     _0_0,      -0.18),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF91_,   epoch=_E(1988),
                              xform=_P(   26.7,    27.5,   -19.9,    2.15,     _0_0,     _0_0,      -0.18),
                              rates=_P(   _0_0,    -0.6,    -1.4,    0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF90_,   epoch=_E(1988),
                              xform=_P(   2.47,     2.35,   -3.59,   2.45,     _0_0,     _0_0,      -0.18),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF89_,   epoch=_E(1988),
                              xform=_P(   2.97,     4.75,   -7.39,   5.85,     _0_0,     _0_0,      -0.18),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))
_trfX(_ITRF2000_, _ITRF88_,   epoch=_E(1988),
                              xform=_P(   2.47,     1.15,   -9.79,   8.95,      0.1,     _0_0,      -0.18),
                              rates=_P(  _0_0,     -0.06,   -0.14,   0.01,     _0_0,     _0_0,       0.02))

# Soler, T .& Snay R.A. U{"Transforming Positions and Velocities between the International Terrestrial
# Reference Frame of 2000 and North American Datum of 1983"<https://www.ResearchGate.net/publication/245291390>},
# Pearson, C. & Snay, R. U{"Introducing HTDP 3.1 to transform coordinates across time and spatial reference
# frames"<https://Geodesy.NOAA.gov/TOOLS/Htdp/Pearson_Snay_2012.pdf> Table 7, 5th and 6th column
_trfX(_ITRF97_,   _ITRF96_,   epoch=_E(1997),
                              xform=_P(  -2.07,    -0.21,    9.95,  -0.93496,   0.1267,   -0.22355,  -0.06065,),
                              rates=_P(   0.69,    -0.1,     1.86,  -0.19201,   0.01347,  -0.01514,   0.00027))
_trfX(_ITRF96_,   _NAD83_,    epoch=_E(1997),
                              xform=_P( 991.0,   -190.72, -512.9,   _0_0,      25.79,      9.65,     11.66,),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.0532,   -0.7423,   -0.0316,))

# see Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG.IFN.Fr/pub/EUREF-TN-1-Jan-31-2024.pdf>} Table 1.
# _trfX(_ITRF2020_, _ETRF2020_, epoch=_E(1989),  # see Table 2 below
#                             xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,      _0_0),
#                             rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.086,     0.519,    -0.753))
# _trfX(_ITRF2014_, _ETRF2014_, epoch=_E(1989),  # see Table 3 below
#                             xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,      _0_0),
#                             rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.085,     0.531,    -0.77))
_trfX(_ITRF2005_, _ETRF2005_, epoch=_E(1989),
                              xform=_P(  56.0,     48.0,   -37.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.054,     0.518,    -0.781))
# _trfX(_ITRF2000_, _ETRF2000_, epoch=_E(1989),  # see Table 4 below
#                             xform=_P(  54.0,     51.0,   -48.0,   _0_0,      _0_0,      _0_0,      _0_0),
#                             rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,     0.49,     -0.792))
_trfX(_ITRF97_,   _ETRF97_,   epoch=_E(1989),
                              xform=_P(  41.0,     41.0,   -49.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,       0.5,      -0.65))
_trfX(_ITRF96_,   _ETRF96_,   epoch=_E(1989),
                              xform=_P(  41.0,     41.0,   -49.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,       0.5,      -0.65))
_trfX(_ITRF94_,   _ETRF94_,   epoch=_E(1989),
                              xform=_P(  41.0,     41.0,   -49.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.2,       0.5,      -0.65))
_trfX(_ITRF93_,   _ETRF93_,   epoch=_E(1989),
                              xform=_P(  19.0,     53.0,   -21.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.32,      0.78,     -0.67))
_trfX(_ITRF92_,   _ETRF92_,   epoch=_E(1989),
                              xform=_P(  38.0,     40.0,   -37.0,    0.0,       0.0,       0.0,       0.0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.21,      0.52,     -0.68))
_trfX(_ITRF91_,   _ETRF91_,   epoch=_E(1989),
                              xform=_P(  21.0,     25.0,   -37.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.21,      0.52,     -0.68))
_trfX(_ITRF90_,   _ETRF90_,   epoch=_E(1989),
                              xform=_P(  19.0,     28.0,   -23.0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.11,      0.57,     -0.71))
_trfX(_ITRF89_,   _ETRF89_,   epoch=_E(1989),
                              xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,      _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.11,      0.57,     -0.71))

# see Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG.IFN.Fr/pub/EUREF-TN-1-Jan-31-2024.pdf>} Table 2.
_trfX(_ITRF2020_, _ETRF2020_, epoch=_E(2015),
                              xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,       2.236,    13.494,   -19.578),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.086,     0.519,    -0.753))
_trfX(_ITRF2014_, _ETRF2020_, epoch=_E(2015),
                              xform=_P(   1.4,      0.9,    -1.4,    0.42,      2.236,    13.494,   -19.578),
                              rates=_P(  _0_0,      0.1,    -0.2,   _0_0,       0.086,     0.519,    -0.753))
_trfX(_ITRF2008_, _ETRF2020_, epoch=_E(2015),
                              xform=_P(   3.0,      2.8,     0.5,    0.55,      2.236,    13.494,   -19.578),
                              rates=_P(  _0_0,      0.1,    -0.3,    0.03,      0.086,     0.519,    -0.753))
_trfX(_ITRF2005_, _ETRF2020_, epoch=_E(2015),
                              xform=_P(   5.5,      1.9,    -4.2,    1.49,      2.236,    13.494,   -19.578),
                              rates=_P(   0.3,      0.1,    -0.3,    0.03,      0.086,     0.519,    -0.753))
_trfX(_ITRF2000_, _ETRF2020_, epoch=_E(2015),
                              xform=_P(   2.6,      2.6,   -37.0,    3.09,      2.236,    13.494,   -19.578),
                              rates=_P(   0.1,      0.2,    -2.1,    0.11,      0.086,     0.519,    -0.753))
_trfX(_ITRF97_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(   9.3,     -2.1,   -80.7,    4.82,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF96_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(   9.3,     -2.1,   -80.7,    4.82,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF94_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(   9.3,     -2.1,   -80.7,    4.82,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF93_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P( -63.0,      3.7,   -74.1,    5.31,     -1.124,     9.164,   -18.828),
                              rates=_P(  -2.8,     _0_0,    -2.7,    0.12,     -0.024,     0.329,    -0.683))
_trfX(_ITRF92_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(  17.3,     -0.1,   -88.7,    4.11,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF91_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(  29.3,     13.9,   -94.7,    5.51,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF90_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(  27.3,      9.9,  -110.7,    5.81,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF89_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(  32.3,     33.9,  -148.7,    9.21,      2.236,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))
_trfX(_ITRF88_,   _ETRF2020_, epoch=_E(2015),
                              xform=_P(  27.3,     -2.1,  -172.7,   12.31,      2.336,    13.494,   -19.218),
                              rates=_P(   0.1,     -0.4,    -3.5,    0.12,      0.086,     0.519,    -0.733))

# see Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG.IFN.Fr/pub/EUREF-TN-1-Jan-31-2024.pdf>} Table 3.
_trfX(_ITRF2020_, _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -1.4,     -0.9,     1.4,    0.42,      2.21,     13.806,   -20.02),
                              rates=_P(  _0_0,     -0.1,     0.2,   _0_0,       0.085,     0.531,    -0.77))
_trfX(_ITRF2014_, _ETRF2014_, epoch=_E(2015),
                              xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,       2.21,     13.806,   -20.02),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.085,     0.531,    -0.77))
_trfX(_ITRF2008_, _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -1.6,     -1.9,    -1.9,   -0.13,      2.21,     13.806,   -20.02),
                              rates=_P(  _0_0,     _0_0,     0.1,   -0.03,      0.085,     0.531,    -0.77))
_trfX(_ITRF2005_, _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -4.1,     -1.0,     2.8,   -1.07,      2.21,     13.806,   -20.02),
                              rates=_P(  -0.3,     _0_0,     0.1,   -0.03,      0.085,     0.531,    -0.77))
_trfX(_ITRF2000_, _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -1.2,     -1.7,    35.6,   -2.67,      2.21,     13.806,   -20.02),
                              rates=_P(  -0.1,     -0.1,     1.9,   -0.11,      0.085,     0.531,    -0.77))
_trfX(_ITRF97_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -7.9,      3.0,    79.3,   -4.40,      2.210,    13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF96_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -7.9,      3.0,    79.3,   -4.40,      2.210,    13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF94_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P(  -7.9,      3.0,    79.3,   -4.40,      2.210,    13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF93_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P(  64.4,     -2.8,    72.7,   -4.89,      5.570,    18.136,   -20.77),
                              rates=_P(   2.8,      0.1,     2.5,   -0.12,      0.195,     0.721,    -0.84))
_trfX(_ITRF92_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P( -15.9,      1.0,    87.3,   -3.69,      2.21,     13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF91_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P( -27.9,    -13.0,    93.3,   -5.09,      2.21,     13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF90_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P( -25.9,     -9.0,   109.3,   -5.39,      2.21,     13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF89_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P( -30.9,    -33.0,   147.3,   -8.79,      2.21,     13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))
_trfX(_ITRF88_,   _ETRF2014_, epoch=_E(2015),
                              xform=_P( -25.9,      3.0,   171.3,  -11.89,      2.11,     13.806,   -20.38),
                              rates=_P(  -0.1,      0.5,     3.3,   -0.12,      0.085,     0.531,    -0.79))

# see U{Altamimi, Z. "EUREF Technical Note 1: Relationship and Transformation between the International and
# the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1-Jan-31-2024.pdf>} Table 4,
# U{Boucher, C. & Altamimi, Z. "Memo: Specifications for reference frame fixing in the analysis of a EUREF GPS
# campaign" (2011) <https://ETRS89.ENSG.IGN.Fr/memo-V8.pdf>} and U{Altamimi, Z. "Key results of ITRF2014 and
# implication to ETRS89 realization", EUREF2016<https://www.EUREF.EU/symposia/2016SanSebastian/01-02-Altamimi.pdf>}.
_trfX(_ITRF2020_, _ETRF2000_, epoch=_E(2015),
                              xform=_P(  53.8,     51.8,   -82.2,    2.25,      2.106,   12.74,    -20.592),
                              rates=_P(   0.1,     _0_0,    -1.7,    0.11,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2014_, _ETRF2000_, epoch=_E(2000),
#                             xform=_P(  53.7,     51.2,   -55.1,    1.02,      0.891,    5.39,     -8.712),
#                             rates=_P(   0.1,      0.1,    -1.9,    0.11,      0.081,    0.49,     -0.792))
_trfX(_ITRF2014_, _ETRF2000_, epoch=_E(2015),
                              xform=_P(  55.2,     52.7,   -83.6,    2.67,      2.106,   12.74,    -20.592),
                              rates=_P(   0.1,      0.1,    -1.9,    0.11,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2008_, _ETRF2000_, epoch=_E(2000),
#                             xform=_P(  52.1,     49.3,   -58.5,    1.34,      0.891,    5.39,     -8.712),
#                             rates=_P(   0.1,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2008_, _ETRF2000_, epoch=_E(2015),
                              xform=_P(  53.6,     50.8,   -85.5,    2.54,      2.106,   12.74,    -20.592),
                              rates=_P(   0.1,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2005_, _ETRF2000_, epoch=_E(2000),
#                             xform=_P(  54.1,     50.2,   -53.8,    0.4,       0.891,    5.39,     -8.712),
#                             rates=_P(  -0.2,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
_trfX(_ITRF2005_, _ETRF2000_, epoch=_E(2015),
                              xform=_P(  51.1,     51.7,   -80.8,    1.6,       2.106,   12.74,    -20.592),
                              rates=_P(  -0.2,      0.1,    -1.8,    0.08,      0.081,    0.49,     -0.792))
# _trfX(_ITRF2000_, _ETRF2000_, epoch=_E(2000),
#                             xform=_P(  54.0,     51.0,   -48.0,   _0_0,       0.891,    5.39,     -8.712),
#                             rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.792))
_trfX(_ITRF2000_, _ETRF2000_, epoch=_E(2015),
                              xform=_P(  54.0,     51.0,   -48.0,   _0_0,       2.106,   12.74,    -20.592),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       0.081,    0.49,     -0.792))
_trfX(_ITRF97_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  47.3,     55.7,    -4.3,   -1.73,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF96_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  47.3,     55.7,    -4.3,   -1.73,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF94_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  47.3,     55.7,    -4.3,   -1.73,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF93_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P( 119.6,     49.9,   -10.9,   -2.22,      5.466,   17.07,    -21.342),
                              rates=_P(   2.9,      0.2,     0.6,   -0.01,      0.191,    0.68,     -0.862))
_trfX(_ITRF92_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  39.3,     53.7,     3.7,   -1.02,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF91_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  27.3,     39.7,     9.7,   -2.42,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF90_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  29.3,     43.7,    25.7,   -2.72,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF89_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  24.3,     19.7,    63.7,   -6.12,      2.106,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))
_trfX(_ITRF88_,   _ETRF2000_, epoch=_E(2015),
                              xform=_P(  29.3,     55.7,    87.7,   -9.22,      2.006,   12.74,    -20.952),
                              rates=_P(  _0_0,      0.6,     1.4,   -0.01,      0.081,    0.49,     -0.812))

# GDA2020 "Geocentric Datum of Australia 2020 Technical Manual", v1.5, 2020-12-09, Table 3.3 and 3.4
# <https://www.ICSM.gov.AU/sites/default/files/2020-12/GDA2020%20Technical%20Manual%20V1.5_4.pdf>
# (the GDA2020 xforms are different but the rates are the same as GDA94, further below)
_trfX(_ITRF2014_, _GDA2020_,  epoch=_E(2020),
                              xform=_P(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,       1.50379,  1.18346,   1.20716))
_trfX(_ITRF2008_, _GDA2020_,  epoch=_E(2020),
                              xform=_P(  13.79,     4.55,   15.22,   2.5,       0.2808,   0.2677,   -0.4638),
                              rates=_P(   1.42,     1.34,    0.9,    0.109,     1.5461,   1.182,     1.1551))
_trfX(_ITRF2005_, _GDA2020_,  epoch=_E(2020),
                              xform=_P(  40.32,   -33.85,  -16.72,   4.286,    -1.2893,  -0.8492,   -0.3342),
                              rates=_P(   2.25,    -0.62,   -0.56,   0.294,    -1.4707,  -1.1443,   -1.1701))
_trfX(_ITRF2000_, _GDA2020_,  epoch=_E(2020),
                              xform=_P(-105.52,    51.58,  231.68,   3.55,      4.2175,   6.3941,    0.8617),
                              rates=_P(  -4.66,     3.55,   11.24,   0.249,     1.7454,   1.4868,    1.224))

# see Table 2 in U{Dawson, J. & Woods, A. "ITRF to GDA94 coordinate transformations", Journal of Applied
# Geodesy 4 (2010), 189-199<https://www.ResearchGate.net/publication/258401581_ITRF_to_GDA94_coordinate_transformations>}
# (note, sign of rotations for GDA94 reversed as "Australia assumes rotation to be of coordinate axes",
# rather than the more conventional "position around the coordinate axes")
_trfX(_ITRF2008_, _GDA94_,    epoch=_E(1994),
                              xform=_P( -84.68,   -19.42,   32.01,   9.71,     -0.4254,   2.2578,    2.4015),
                              rates=_P(   1.42,     1.34,    0.9,    0.109,     1.5461,   1.182,     1.1551))
_trfX(_ITRF2005_, _GDA94_,    epoch=_E(1994),
                              xform=_P( -79.73,    -6.86,   38.03,   6.636,     0.0351,  -2.1211,   -2.1411),
                              rates=_P(   2.25,    -0.62,   -0.56,   0.294,    -1.4707,  -1.1443,   -1.1701))
_trfX(_ITRF2000_, _GDA94_,    epoch=_E(1994),
                              xform=_P( -45.91,   -29.85,  -20.37,   7.07,     -1.6705,   0.4594,    1.9356),
                              rates=_P(  -4.66,     3.55,   11.24,   0.249,     1.7454,   1.4868,    1.224))

# see U{Solar, T. & Snay, R.A. "Transforming Positions and Velocities between the
# International Terrestrial Reference Frame of 2000 and North American Datum of 1983"
# (2004)<https://www.NGS.NOAA.gov/CORS/Articles/SolerSnayASCE.pdf>}
_trfX(_ITRF2000_, _NAD83_,    epoch=_E(1997),  # note NAD83(CORS96)
                              xform=_P( 995.6,  -1901.3,  -521.5,    0.615,   -25.915,   -9.426,   -11.599),
                              rates=_P(   0.7,     -0.7,    _0_5,   -0.182,    -0.06667,  0.75744,   0.05133))

# see U{Quinsy QPS<https://confluence.QPS.NL/qinsy/files/latest/en/182618383/182618384/1/1579182881000/
# ITRF_Transformation_Parameters.xlsx>}, sheets ITRF and NAD83 and Pearson, C. & Snay, R. U{"Introducing
# HTDP 3.1 to transform coordinates across time and spatial reference frames"<https://Geodesy.NOAA.gov/
# TOOLS/Htdp/Pearson_Snay_2012.pdf> Table 7, 7th column
_trfX(_ITRF2008_, _NAD83_,    epoch=_E(1997),
                              xform=_P( 993.43, -1903.31, -526.55,   1.71504, -25.91467, -9.42645, -11.59935),
                              rates=_P(   0.79,    -0.6,    -1.34,  -0.10201,  -0.06667,  0.75744,   0.05133))
# see U{Quinsy QPS<https://confluence.QPS.NL/qinsy/files/latest/en/182618383/182618384/1/1579182881000/
# ITRF_Transformation_Parameters.xlsx>}, sheets ITRF and NAD83
_trfX(_ITRF2005_, _NAD83_,    epoch=_E(1997),
                              xform=_P( 996.3,  -1902.4,  -521.9,    0.775,   -25.915,   -9.426,   -11.599),
                              rates=_P(   0.5,     -0.6,    -1.3,   -0.10201,  -0.06667,  0.75744,   0.05133))
_trfX(_ITRF90_,   _NAD83_,    epoch=_E(1997),
                              xform=_P( 973.0,  -1919.2,  -482.9,   -0.9,     -25.79,    -9.65,    -11.66),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,      -0.053,    0.742,     0.032))
_trfX(_ITRF90_,   _WGS84_,    epoch=_E(1984),
                              xform=_P(  60.0,   -517.0,  -223.0,  -11.0,      18.3,     -0.3,       7.0),
                              rates=_P(  _0_0,     _0_0,    _0_0,   _0_0,      _0_0,     _0_0,      _0_0))

del _P, _Ps  # keep _E, _Es

if __name__ == '__main__':

    from pygeodesy.interns import _COLONSPACE_,_COMMA_, _NL_, _NLATvar_, _STAR_
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
    def _RFs():
        yield NN
        for t in RefFrames.toRepr(all=True).split(_NL_):
            t = t.strip(_COMMA_)
            n = t.split(_COLONSPACE_)[0].split(_DOT_)[1]
            r = RefFrames.get(n)
            x = len(r.Xforms()), -len(r.Xforms(inverse=True))
            yield '%s .Xforms=%s' % (t, x)

    printf(_NLATvar_.join(sorted(_RFs())), nt=1)

    X, t = (), 0  # all  form
    for r in RefFrames.values():
        X += tuple(r._Xs.values())
    for X in sorted(X):
        t += 1
        printf('#%4d %-27s xform=%r', t, X.name, X.xform)
        printf('#%32s rates=%r', _SPACE_, X.rates)

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

#   1 ITRF2000@1988xITRF88@1988   xform=(tx=2.47, ty=1.15, tz=-9.79, s=8.95, sx=0.1, sy=0.0, sz=-0.18)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   2 ITRF2000@1988xITRF89@1988   xform=(tx=2.97, ty=4.75, tz=-7.39, s=5.85, sx=0.0, sy=0.0, sz=-0.18)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   3 ITRF2000@1988xITRF90@1988   xform=(tx=2.47, ty=2.35, tz=-3.59, s=2.45, sx=0.0, sy=0.0, sz=-0.18)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   4 ITRF2000@1988xITRF91@1988   xform=(tx=26.7, ty=27.5, tz=-19.9, s=2.15, sx=0.0, sy=0.0, sz=-0.18)
#                                 rates=(tx=0.0, ty=-0.6, tz=-1.4, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   5 ITRF2000@1988xITRF92@1988   xform=(tx=1.47, ty=1.35, tz=-1.39, s=0.75, sx=0.0, sy=0.0, sz=-0.18)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   6 ITRF2000@1988xITRF93@1988   xform=(tx=12.7, ty=6.5, tz=-20.9, s=1.95, sx=-0.39, sy=0.8, sz=-1.14)
#                                 rates=(tx=-2.9, ty=-0.2, tz=-0.6, s=0.01, sx=-0.11, sy=-0.19, sz=0.07)
#   7 ITRF2000@1994xGDA94@1994    xform=(tx=-45.91, ty=-29.85, tz=-20.37, s=7.07, sx=-1.6705, sy=0.4594, sz=1.9356)
#                                 rates=(tx=-4.66, ty=3.55, tz=11.24, s=0.249, sx=1.7454, sy=1.4868, sz=1.224)
#   8 ITRF2000@1997xITRF94@1997   xform=(tx=0.67, ty=0.61, tz=-1.85, s=1.55, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#   9 ITRF2000@1997xITRF96@1997   xform=(tx=0.67, ty=0.61, tz=-1.85, s=1.55, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=0.02)
#  10 ITRF2000@1997xITRF97@1997   xform=(tx=0.67, ty=0.61, tz=-1.85, s=1.55, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=-0.06, tz=-0.14, s=0.01, sx=0.0, sy=0.0, sz=-0.02)
#  11 ITRF2000@1997xNAD83@1997    xform=(tx=995.6, ty=-1901.3, tz=-521.5, s=0.615, sx=-25.915, sy=-9.426, sz=-11.599)
#                                 rates=(tx=0.7, ty=-0.7, tz=0.5, s=-0.182, sx=-0.06667, sy=0.75744, sz=0.05133)
#  12 ITRF2000@2015xETRF2000@2015 xform=(tx=54.0, ty=51.0, tz=-48.0, s=0.0, sx=2.106, sy=12.74, sz=-20.592)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.081, sy=0.49, sz=-0.792)
#  13 ITRF2000@2015xETRF2014@2015 xform=(tx=-1.2, ty=-1.7, tz=35.6, s=-2.67, sx=2.21, sy=13.806, sz=-20.02)
#                                 rates=(tx=-0.1, ty=-0.1, tz=1.9, s=-0.11, sx=0.085, sy=0.531, sz=-0.77)
#  14 ITRF2000@2015xETRF2020@2015 xform=(tx=2.6, ty=2.6, tz=-37.0, s=3.09, sx=2.236, sy=13.494, sz=-19.578)
#                                 rates=(tx=0.1, ty=0.2, tz=-2.1, s=0.11, sx=0.086, sy=0.519, sz=-0.753)
#  15 ITRF2000@2020xGDA2020@2020  xform=(tx=-105.52, ty=51.58, tz=231.68, s=3.55, sx=4.2175, sy=6.3941, sz=0.8617)
#                                 rates=(tx=-4.66, ty=3.55, tz=11.24, s=0.249, sx=1.7454, sy=1.4868, sz=1.224)
#  16 ITRF2005@1989xETRF2005@1989 xform=(tx=56.0, ty=48.0, tz=-37.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.054, sy=0.518, sz=-0.781)
#  17 ITRF2005@1994xGDA94@1994    xform=(tx=-79.73, ty=-6.86, tz=38.03, s=6.636, sx=0.0351, sy=-2.1211, sz=-2.1411)
#                                 rates=(tx=2.25, ty=-0.62, tz=-0.56, s=0.294, sx=-1.4707, sy=-1.1443, sz=-1.1701)
#  18 ITRF2005@1997xNAD83@1997    xform=(tx=996.3, ty=-1902.4, tz=-521.9, s=0.775, sx=-25.915, sy=-9.426, sz=-11.599)
#                                 rates=(tx=0.5, ty=-0.6, tz=-1.3, s=-0.10201, sx=-0.06667, sy=0.75744, sz=0.05133)
#  19 ITRF2005@2000xITRF2000@2000 xform=(tx=0.1, ty=-0.8, tz=-5.8, s=0.4, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=-0.2, ty=0.1, tz=-1.8, s=0.08, sx=0.0, sy=0.0, sz=0.0)
#  20 ITRF2005@2015xETRF2000@2015 xform=(tx=51.1, ty=51.7, tz=-80.8, s=1.6, sx=2.106, sy=12.74, sz=-20.592)
#                                 rates=(tx=-0.2, ty=0.1, tz=-1.8, s=0.08, sx=0.081, sy=0.49, sz=-0.792)
#  21 ITRF2005@2015xETRF2014@2015 xform=(tx=-4.1, ty=-1.0, tz=2.8, s=-1.07, sx=2.21, sy=13.806, sz=-20.02)
#                                 rates=(tx=-0.3, ty=0.0, tz=0.1, s=-0.03, sx=0.085, sy=0.531, sz=-0.77)
#  22 ITRF2005@2015xETRF2020@2015 xform=(tx=5.5, ty=1.9, tz=-4.2, s=1.49, sx=2.236, sy=13.494, sz=-19.578)
#                                 rates=(tx=0.3, ty=0.1, tz=-0.3, s=0.03, sx=0.086, sy=0.519, sz=-0.753)
#  23 ITRF2005@2020xGDA2020@2020  xform=(tx=40.32, ty=-33.85, tz=-16.72, s=4.286, sx=-1.2893, sy=-0.8492, sz=-0.3342)
#                                 rates=(tx=2.25, ty=-0.62, tz=-0.56, s=0.294, sx=-1.4707, sy=-1.1443, sz=-1.1701)
#  24 ITRF2008@1994xGDA94@1994    xform=(tx=-84.68, ty=-19.42, tz=32.01, s=9.71, sx=-0.4254, sy=2.2578, sz=2.4015)
#                                 rates=(tx=1.42, ty=1.34, tz=0.9, s=0.109, sx=1.5461, sy=1.182, sz=1.1551)
#  25 ITRF2008@1997xNAD83@1997    xform=(tx=993.43, ty=-1903.31, tz=-526.55, s=1.71504, sx=-25.91467, sy=-9.42645, sz=-11.59935)
#                                 rates=(tx=0.79, ty=-0.6, tz=-1.34, s=-0.10201, sx=-0.06667, sy=0.75744, sz=0.05133)
#  26 ITRF2008@2000xITRF2000@2000 xform=(tx=-1.9, ty=-1.7, tz=-10.5, s=1.34, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.1, ty=0.1, tz=-1.8, s=0.08, sx=0.0, sy=0.0, sz=0.0)
#  27 ITRF2008@2000xITRF88@2000   xform=(tx=22.8, ty=2.6, tz=-125.2, s=10.41, sx=0.1, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  28 ITRF2008@2000xITRF89@2000   xform=(tx=27.8, ty=38.6, tz=-101.2, s=7.31, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  29 ITRF2008@2000xITRF90@2000   xform=(tx=22.8, ty=14.6, tz=-63.2, s=3.91, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  30 ITRF2008@2000xITRF91@2000   xform=(tx=24.8, ty=18.6, tz=-47.2, s=3.61, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  31 ITRF2008@2000xITRF92@2000   xform=(tx=12.8, ty=4.6, tz=-41.2, s=2.21, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  32 ITRF2008@2000xITRF93@2000   xform=(tx=-24.0, ty=2.4, tz=-38.6, s=3.41, sx=-1.71, sy=-1.48, sz=-0.3)
#                                 rates=(tx=-2.8, ty=-0.1, tz=-2.4, s=0.09, sx=-0.11, sy=-0.19, sz=0.07)
#  33 ITRF2008@2000xITRF94@2000   xform=(tx=4.8, ty=2.6, tz=-33.2, s=2.92, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  34 ITRF2008@2000xITRF96@2000   xform=(tx=4.8, ty=2.6, tz=-33.2, s=2.92, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  35 ITRF2008@2000xITRF97@2000   xform=(tx=4.8, ty=2.6, tz=-33.2, s=2.92, sx=0.0, sy=0.0, sz=0.06)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.2, s=0.09, sx=0.0, sy=0.0, sz=0.02)
#  36 ITRF2008@2005xITRF2005@2005 xform=(tx=-0.5, ty=-0.9, tz=-4.7, s=0.94, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.3, ty=0.0, tz=0.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#  37 ITRF2008@2015xETRF2000@2015 xform=(tx=53.6, ty=50.8, tz=-85.5, s=2.54, sx=2.106, sy=12.74, sz=-20.592)
#                                 rates=(tx=0.1, ty=0.1, tz=-1.8, s=0.08, sx=0.081, sy=0.49, sz=-0.792)
#  38 ITRF2008@2015xETRF2014@2015 xform=(tx=-1.6, ty=-1.9, tz=-1.9, s=-0.13, sx=2.21, sy=13.806, sz=-20.02)
#                                 rates=(tx=0.0, ty=0.0, tz=0.1, s=-0.03, sx=0.085, sy=0.531, sz=-0.77)
#  39 ITRF2008@2015xETRF2020@2015 xform=(tx=3.0, ty=2.8, tz=0.5, s=0.55, sx=2.236, sy=13.494, sz=-19.578)
#                                 rates=(tx=0.0, ty=0.1, tz=-0.3, s=0.03, sx=0.086, sy=0.519, sz=-0.753)
#  40 ITRF2008@2020xGDA2020@2020  xform=(tx=13.79, ty=4.55, tz=15.22, s=2.5, sx=0.2808, sy=0.2677, sz=-0.4638)
#                                 rates=(tx=1.42, ty=1.34, tz=0.9, s=0.109, sx=1.5461, sy=1.182, sz=1.1551)
#  41 ITRF2014@2010xITRF2000@2010 xform=(tx=0.7, ty=1.2, tz=-26.1, s=2.12, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.1, ty=0.1, tz=-1.9, s=0.11, sx=0.0, sy=0.0, sz=0.0)
#  42 ITRF2014@2010xITRF2005@2010 xform=(tx=2.6, ty=1.0, tz=-2.3, s=0.92, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.3, ty=0.0, tz=-0.1, s=0.03, sx=0.0, sy=0.0, sz=0.0)
#  43 ITRF2014@2010xITRF2008@2010 xform=(tx=1.6, ty=1.9, tz=2.4, s=-0.02, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=-0.1, s=0.03, sx=0.0, sy=0.0, sz=0.0)
#  44 ITRF2014@2010xITRF88@2010   xform=(tx=25.4, ty=-0.5, tz=-154.8, s=11.29, sx=0.1, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  45 ITRF2014@2010xITRF89@2010   xform=(tx=30.4, ty=35.5, tz=-130.8, s=8.19, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  46 ITRF2014@2010xITRF90@2010   xform=(tx=25.4, ty=11.5, tz=-92.8, s=4.79, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  47 ITRF2014@2010xITRF91@2010   xform=(tx=27.4, ty=15.5, tz=-76.8, s=4.49, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  48 ITRF2014@2010xITRF92@2010   xform=(tx=15.4, ty=1.5, tz=-70.8, s=3.09, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  49 ITRF2014@2010xITRF93@2010   xform=(tx=-50.4, ty=3.3, tz=-60.2, s=4.29, sx=-2.81, sy=-3.38, sz=0.4)
#                                 rates=(tx=-2.8, ty=-0.1, tz=-2.5, s=0.12, sx=-0.11, sy=-0.19, sz=0.07)
#  50 ITRF2014@2010xITRF94@2010   xform=(tx=7.4, ty=-0.5, tz=-62.8, s=3.8, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  51 ITRF2014@2010xITRF96@2010   xform=(tx=7.4, ty=-0.5, tz=-62.8, s=3.8, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  52 ITRF2014@2010xITRF97@2010   xform=(tx=7.4, ty=-0.5, tz=-62.8, s=3.8, sx=0.0, sy=0.0, sz=0.26)
#                                 rates=(tx=0.1, ty=-0.5, tz=-3.3, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  53 ITRF2014@2015xETRF2000@2015 xform=(tx=55.2, ty=52.7, tz=-83.6, s=2.67, sx=2.106, sy=12.74, sz=-20.592)
#                                 rates=(tx=0.1, ty=0.1, tz=-1.9, s=0.11, sx=0.081, sy=0.49, sz=-0.792)
#  54 ITRF2014@2015xETRF2014@2015 xform=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=2.21, sy=13.806, sz=-20.02)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.085, sy=0.531, sz=-0.77)
#  55 ITRF2014@2015xETRF2020@2015 xform=(tx=1.4, ty=0.9, tz=-1.4, s=0.42, sx=2.236, sy=13.494, sz=-19.578)
#                                 rates=(tx=0.0, ty=0.1, tz=-0.2, s=0.0, sx=0.086, sy=0.519, sz=-0.753)
#  56 ITRF2014@2020xGDA2020@2020  xform=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=1.50379, sy=1.18346, sz=1.20716)
#  57 ITRF2020@2015xETRF2000@2015 xform=(tx=53.8, ty=51.8, tz=-82.2, s=2.25, sx=2.106, sy=12.74, sz=-20.592)
#                                 rates=(tx=0.1, ty=0.0, tz=-1.7, s=0.11, sx=0.081, sy=0.49, sz=-0.792)
#  58 ITRF2020@2015xETRF2014@2015 xform=(tx=-1.4, ty=-0.9, tz=1.4, s=0.42, sx=2.21, sy=13.806, sz=-20.02)
#                                 rates=(tx=0.0, ty=-0.1, tz=0.2, s=0.0, sx=0.085, sy=0.531, sz=-0.77)
#  59 ITRF2020@2015xETRF2020@2015 xform=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=2.236, sy=13.494, sz=-19.578)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.086, sy=0.519, sz=-0.753)
#  60 ITRF2020@2015xITRF2000@2015 xform=(tx=-0.2, ty=0.8, tz=-34.2, s=2.25, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.1, ty=0.0, tz=-1.7, s=0.11, sx=0.0, sy=0.0, sz=0.0)
#  61 ITRF2020@2015xITRF2005@2015 xform=(tx=2.7, ty=0.1, tz=-1.4, s=0.65, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.3, ty=-0.1, tz=0.1, s=0.03, sx=0.0, sy=0.0, sz=0.0)
#  62 ITRF2020@2015xITRF2008@2015 xform=(tx=0.2, ty=1.0, tz=3.3, s=-0.29, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=-0.1, tz=0.1, s=0.03, sx=0.0, sy=0.0, sz=0.0)
#  63 ITRF2020@2015xITRF2014@2015 xform=(tx=-1.4, ty=-0.9, tz=1.4, s=-0.42, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=-0.1, tz=0.2, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#  64 ITRF2020@2015xITRF88@2015   xform=(tx=24.5, ty=-3.9, tz=-169.9, s=11.47, sx=0.1, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  65 ITRF2020@2015xITRF89@2015   xform=(tx=29.5, ty=32.1, tz=-145.9, s=8.37, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  66 ITRF2020@2015xITRF90@2015   xform=(tx=24.5, ty=8.1, tz=-107.9, s=4.97, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  67 ITRF2020@2015xITRF91@2015   xform=(tx=26.5, ty=12.1, tz=-91.9, s=4.67, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  68 ITRF2020@2015xITRF92@2015   xform=(tx=14.5, ty=-1.9, tz=-85.9, s=3.27, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  69 ITRF2020@2015xITRF93@2015   xform=(tx=-65.8, ty=1.9, tz=-71.3, s=4.47, sx=-3.36, sy=-4.33, sz=0.75)
#                                 rates=(tx=-2.8, ty=-0.2, tz=-2.3, s=0.12, sx=-0.11, sy=-0.19, sz=0.07)
#  70 ITRF2020@2015xITRF94@2015   xform=(tx=6.5, ty=-3.9, tz=-77.9, s=3.98, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  71 ITRF2020@2015xITRF96@2015   xform=(tx=6.5, ty=-3.9, tz=-77.9, s=3.98, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  72 ITRF2020@2015xITRF97@2015   xform=(tx=6.5, ty=-3.9, tz=-77.9, s=3.98, sx=0.0, sy=0.0, sz=0.36)
#                                 rates=(tx=0.1, ty=-0.6, tz=-3.1, s=0.12, sx=0.0, sy=0.0, sz=0.02)
#  73 ITRF88@2015xETRF2000@2015   xform=(tx=29.3, ty=55.7, tz=87.7, s=-9.22, sx=2.006, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
#  74 ITRF88@2015xETRF2014@2015   xform=(tx=-25.9, ty=3.0, tz=171.3, s=-11.89, sx=2.11, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
#  75 ITRF88@2015xETRF2020@2015   xform=(tx=27.3, ty=-2.1, tz=-172.7, s=12.31, sx=2.336, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
#  76 ITRF89@1989xETRF89@1989     xform=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.11, sy=0.57, sz=-0.71)
#  77 ITRF89@2015xETRF2000@2015   xform=(tx=24.3, ty=19.7, tz=63.7, s=-6.12, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
#  78 ITRF89@2015xETRF2014@2015   xform=(tx=-30.9, ty=-33.0, tz=147.3, s=-8.79, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
#  79 ITRF89@2015xETRF2020@2015   xform=(tx=32.3, ty=33.9, tz=-148.7, s=9.21, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
#  80 ITRF90@1984xWGS84@1984      xform=(tx=60.0, ty=-517.0, tz=-223.0, s=-11.0, sx=18.3, sy=-0.3, sz=7.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#  81 ITRF90@1989xETRF90@1989     xform=(tx=19.0, ty=28.0, tz=-23.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.11, sy=0.57, sz=-0.71)
#  82 ITRF90@1997xNAD83@1997      xform=(tx=973.0, ty=-1919.2, tz=-482.9, s=-0.9, sx=-25.79, sy=-9.65, sz=-11.66)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=-0.053, sy=0.742, sz=0.032)
#  83 ITRF90@2015xETRF2000@2015   xform=(tx=29.3, ty=43.7, tz=25.7, s=-2.72, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
#  84 ITRF90@2015xETRF2014@2015   xform=(tx=-25.9, ty=-9.0, tz=109.3, s=-5.39, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
#  85 ITRF90@2015xETRF2020@2015   xform=(tx=27.3, ty=9.9, tz=-110.7, s=5.81, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
#  86 ITRF91@1989xETRF91@1989     xform=(tx=21.0, ty=25.0, tz=-37.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.21, sy=0.52, sz=-0.68)
#  87 ITRF91@2015xETRF2000@2015   xform=(tx=27.3, ty=39.7, tz=9.7, s=-2.42, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
#  88 ITRF91@2015xETRF2014@2015   xform=(tx=-27.9, ty=-13.0, tz=93.3, s=-5.09, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
#  89 ITRF91@2015xETRF2020@2015   xform=(tx=29.3, ty=13.9, tz=-94.7, s=5.51, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
#  90 ITRF92@1989xETRF92@1989     xform=(tx=38.0, ty=40.0, tz=-37.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.21, sy=0.52, sz=-0.68)
#  91 ITRF92@2015xETRF2000@2015   xform=(tx=39.3, ty=53.7, tz=3.7, s=-1.02, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
#  92 ITRF92@2015xETRF2014@2015   xform=(tx=-15.9, ty=1.0, tz=87.3, s=-3.69, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
#  93 ITRF92@2015xETRF2020@2015   xform=(tx=17.3, ty=-0.1, tz=-88.7, s=4.11, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
#  94 ITRF93@1989xETRF93@1989     xform=(tx=19.0, ty=53.0, tz=-21.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.32, sy=0.78, sz=-0.67)
#  95 ITRF93@2015xETRF2000@2015   xform=(tx=119.6, ty=49.9, tz=-10.9, s=-2.22, sx=5.466, sy=17.07, sz=-21.342)
#                                 rates=(tx=2.9, ty=0.2, tz=0.6, s=-0.01, sx=0.191, sy=0.68, sz=-0.862)
#  96 ITRF93@2015xETRF2014@2015   xform=(tx=64.4, ty=-2.8, tz=72.7, s=-4.89, sx=5.57, sy=18.136, sz=-20.77)
#                                 rates=(tx=2.8, ty=0.1, tz=2.5, s=-0.12, sx=0.195, sy=0.721, sz=-0.84)
#  97 ITRF93@2015xETRF2020@2015   xform=(tx=-63.0, ty=3.7, tz=-74.1, s=5.31, sx=-1.124, sy=9.164, sz=-18.828)
#                                 rates=(tx=-2.8, ty=0.0, tz=-2.7, s=0.12, sx=-0.024, sy=0.329, sz=-0.683)
#  98 ITRF94@1989xETRF94@1989     xform=(tx=41.0, ty=41.0, tz=-49.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.2, sy=0.5, sz=-0.65)
#  99 ITRF94@2015xETRF2000@2015   xform=(tx=47.3, ty=55.7, tz=-4.3, s=-1.73, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
# 100 ITRF94@2015xETRF2014@2015   xform=(tx=-7.9, ty=3.0, tz=79.3, s=-4.4, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
# 101 ITRF94@2015xETRF2020@2015   xform=(tx=9.3, ty=-2.1, tz=-80.7, s=4.82, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
# 102 ITRF96@1989xETRF96@1989     xform=(tx=41.0, ty=41.0, tz=-49.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.2, sy=0.5, sz=-0.65)
# 103 ITRF96@1997xNAD83@1997      xform=(tx=991.0, ty=-190.72, tz=-512.9, s=0.0, sx=25.79, sy=9.65, sz=11.66)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.0532, sy=-0.7423, sz=-0.0316)
# 104 ITRF96@2015xETRF2000@2015   xform=(tx=47.3, ty=55.7, tz=-4.3, s=-1.73, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
# 105 ITRF96@2015xETRF2014@2015   xform=(tx=-7.9, ty=3.0, tz=79.3, s=-4.4, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
# 106 ITRF96@2015xETRF2020@2015   xform=(tx=9.3, ty=-2.1, tz=-80.7, s=4.82, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
# 107 ITRF97@1989xETRF97@1989     xform=(tx=41.0, ty=41.0, tz=-49.0, s=0.0, sx=0.0, sy=0.0, sz=0.0)
#                                 rates=(tx=0.0, ty=0.0, tz=0.0, s=0.0, sx=0.2, sy=0.5, sz=-0.65)
# 108 ITRF97@1997xITRF96@1997     xform=(tx=-2.07, ty=-0.21, tz=9.95, s=-0.93496, sx=0.1267, sy=-0.22355, sz=-0.06065)
#                                 rates=(tx=0.69, ty=-0.1, tz=1.86, s=-0.19201, sx=0.01347, sy=-0.01514, sz=0.00027)
# 109 ITRF97@2015xETRF2000@2015   xform=(tx=47.3, ty=55.7, tz=-4.3, s=-1.73, sx=2.106, sy=12.74, sz=-20.952)
#                                 rates=(tx=0.0, ty=0.6, tz=1.4, s=-0.01, sx=0.081, sy=0.49, sz=-0.812)
# 110 ITRF97@2015xETRF2014@2015   xform=(tx=-7.9, ty=3.0, tz=79.3, s=-4.4, sx=2.21, sy=13.806, sz=-20.38)
#                                 rates=(tx=-0.1, ty=0.5, tz=3.3, s=-0.12, sx=0.085, sy=0.531, sz=-0.79)
# 111 ITRF97@2015xETRF2020@2015   xform=(tx=9.3, ty=-2.1, tz=-80.7, s=4.82, sx=2.236, sy=13.494, sz=-19.218)
#                                 rates=(tx=0.1, ty=-0.4, tz=-3.5, s=0.12, sx=0.086, sy=0.519, sz=-0.733)
