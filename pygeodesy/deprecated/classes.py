
# -*- coding: utf-8 -*-

u'''DEPRECATED classes for export and backward compatibility.
'''

from pygeodesy.clipy import ClipCS4Tuple as _ClipCS4Tuple
from pygeodesy.constants import NAN, _float
from pygeodesy.interns import NN, _a12_, _area_, _band_, _convergence_, \
                             _gamma_, _i_, _ltp_
from pygeodesy.deprecated.consterns import _Deprecated_Str
from pygeodesy.karney import Rhumb8Tuple as _Rhumb8Tuple
from pygeodesy.lazily import _ALL_MODS as _MODS, _ALL_OTHER
from pygeodesy.ltpTuples import Ned4Tuple as _Ned4Tuple
# from pygeodesy.named import _NamedTuple  # from .namedTuples
from pygeodesy.namedTuples import Forward4Tuple as _Forward4Tuple, \
                                  Reverse4Tuple as _Reverse4Tuple, \
                                  UtmUps5Tuple as _UtmUps5Tuple,  _NamedTuple
from pygeodesy.props import deprecated_class, deprecated_method
from pygeodesy.resections import TriAngle5Tuple as _TriAngle5Tuple
from pygeodesy.trf import Helmert7Tuple as _Helmert7Tuple

__all__ = ()
__version__ = '23.09.12'


class _Deprecated_NamedTuple(_NamedTuple):
    '''DEPRECATED, C{_NamedTuple} base.
    '''
    def __new__(cls, *args, **kwds):
        deprecated_class(cls)
        return _NamedTuple.__new__(cls, *args, **kwds)


def _reNames(names, old, *new):
    # replace item C{old} with C{new} name
    i = names.index(old)
    return names[:i] + new + names[i + 1:]


class ClipCS3Tuple(_Deprecated_NamedTuple):  # PYCHOK no cover
    '''DEPRECATED, see I{DEPRECATED} function L{pygeodesy.clipCS3}.'''
    assert _ClipCS4Tuple._Names_.index(_i_) == 2
    _Names_ = _reNames(_ClipCS4Tuple._Names_[:3], _i_, 'index')
    _Units_ =          _ClipCS4Tuple._Units_[:3]


class EasNorExact4Tuple(_Deprecated_NamedTuple):
    '''DEPRECATED, use class L{Forward4Tuple}, item C{gamma} for C{convergence}.'''
    _Names_ = _reNames(_Forward4Tuple._Names_, _gamma_, _convergence_)
    _Units_ =          _Forward4Tuple._Units_


def EcefCartesian(*args, **kwds):
    '''DEPRECATED, use class L{LocalCartesian}.'''
    Ltp = _MODS.ltp.Ltp

    class EcefCartesian_(Ltp):
        '''DEPRECATED, use class L{LocalCartesian} or L{Ltp}.

           @note: This class is named I{incorrectly}, since it provides conversion to
                  and from I{local} cartesian coordinates in a I{local tangent plane}
                  and I{not geocentric} (ECEF) ones, as the name would suggest.
        '''
        def __init__(self, latlonh0=0, lon0=0, height0=0, ecef=None, name=NN):
            deprecated_class(self.__class__)
            Ltp.__init__(self, latlonh0=latlonh0, lon0=lon0, height0=height0, ecef=ecef, name=name)

        @deprecated_method
        def forward(self, latlonh, lon=None, height=0, M=False, name=NN):
            '''DEPRECATED, use method L{LocalCartesian.forward} or L{Ltp.forward}.

               @return: I{Incorrectly}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C,
                        M, datum)} with I{local} C{(x, y, z)} coordinates for the given
                        I{geodetic} ones C{(lat, lon, height)}, case C{C=0} always,
                        optionally I{concatenated} L{EcefMatrix} C{M} and C{datum}.
            '''
            t = Ltp.forward(self, latlonh, lon=lon, height=height, M=M, name=name)
            return _MODS.ecef.Ecef9Tuple(t.x, t.y, t.z, t.lat, t.lon, t.height,
                                                        0, t.M, t.ecef.datum,
                                                        name=t.name or self.name)

        @deprecated_method
        def reverse(self, xyz, y=None, z=None, M=False, name=NN):
            '''DEPRECATED, use method L{LocalCartesian.reverse} or L{Ltp.reverse}.

               @return: I{Incorrectly}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C,
                        M, datum)} with I{geodetic} coordinates C{(lat, lon, height)} for
                        the given I{local} ones C{(x, y, z)}, case C{C}, optionally
                        I{concatenated} L{EcefMatrix} C{M} and C{datum}.
            '''
            t = Ltp.reverse(self, xyz, y=y, z=z, M=M, name=name)
            return _MODS.ecef.Ecef9Tuple(t.x, t.y, t.z, t.lat, t.lon, t.height,
                                                        t.ecef.C, t.M, t.ecef.datum,
                                                        name=t.name or self.name)

    _MODS.deprecated.EcefCartesian = Ltp
    return EcefCartesian_(*args, **kwds)


def HeightIDW(knots, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use class L{HeightIDWeuclidean}.'''
    HeightIDWeuclidean = _MODS.heights.HeightIDWeuclidean

    class HeightIDW(HeightIDWeuclidean):
        '''DEPRECATED, use class L{HeightIDWeuclidean}.'''
        def __init__(self, knots, adjust=True, beta=2, name=NN):
            deprecated_class(self.__class__)
            HeightIDWeuclidean.__init__(self, knots, adjust=adjust, beta=beta, name=name)

    _MODS.deprecated.HeightIDW = HeightIDW
    return HeightIDW(knots, **kwds)


def HeightIDW2(knots, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use class L{HeightIDWequirectangular}.'''
    HeightIDWequirectangular = _MODS.heights.HeightIDWequirectangular

    class HeightIDW2(HeightIDWequirectangular):
        '''DEPRECATED, use class L{HeightIDWequirectangular}.'''
        def __init__(self, knots, adjust=True, wrap=False, name=NN):
            deprecated_class(self.__class__)
            HeightIDWequirectangular.__init__(self, knots, adjust=adjust, wrap=wrap, name=name)

    _MODS.deprecated.HeightIDW2 = HeightIDW2
    return HeightIDW2(knots, **kwds)


def HeightIDW3(knots, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use class L{HeightIDWhaversine}.'''
    HeightIDWhaversine = _MODS.heights.HeightIDWhaversine

    class HeightIDW3(HeightIDWhaversine):
        '''DEPRECATED, use class L{HeightIDWhaversine}.
        '''
        def __init__(self, knots, beta=2, wrap=False, name=NN):
            deprecated_class(self.__class__)
            HeightIDWhaversine.__init__(self, knots, beta=beta, wrap=wrap, name=name)

    _MODS.deprecated.HeightIDW3 = HeightIDW3
    return HeightIDW3(knots, **kwds)


class LatLonExact4Tuple(_Deprecated_NamedTuple):
    '''DEPRECATED, used class L{Reverse4Tuple}, item C{gamma} for C{convergence}.'''
    _Names_ = _reNames(_Reverse4Tuple._Names_, _gamma_, _convergence_)
    _Units_ =          _Reverse4Tuple._Units_


class Ned3Tuple(_Deprecated_NamedTuple):  # was in .ellipsoidalNvector
    '''DEPRECATED, use class L{Ned4Tuple}, ignoring item C{ltp}.'''
    assert _Ned4Tuple._Names_.index(_ltp_) == 3
    _Names_ = _Ned4Tuple._Names_[:3]
    _Units_ = _Ned4Tuple._Units_[:3]


def RefFrameError(*args, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use class L{TRFError}.'''
    TRFError = _MODS.errors.TRFError

    class RefFrameError(TRFError):
        '''DEPRECATED, use class L{TRFError}.
        '''
        def __init__(self, *name_value, **txt_name_values):
            deprecated_class(self.__class__)
            TRFError.__init__(self, *name_value, **txt_name_values)

    _MODS.deprecated.RefFrameError = RefFrameError
    return RefFrameError(*args, **kwds)


class Rhumb7Tuple(_Deprecated_NamedTuple):
    '''DEPRECATED, use class L{Rhumb8Tuple}, ignoring item C{a12}.'''
    assert _Rhumb8Tuple._Names_.index(_a12_) == 7
    _Names_ = _Rhumb8Tuple._Names_[:7]
    _Units_ = _Rhumb8Tuple._Units_[:7]

    @deprecated_method
    def toDirect9Tuple(self, **kwds):
        return self.toRhumb8Tuple().toDirect9Tuple(self, **kwds)

    @deprecated_method
    def toGDict(self, **kwds):
        return self.toRhumb8Tuple().toGDict(**kwds)

    @deprecated_method
    def toInverse10Tuple(self, **kwds):
        return self.toRhumb8Tuple().toInverse10Tuple(self, **kwds)

    @deprecated_method
    def toRhumb8Tuple(self, dflt=NAN):
        return _Rhumb8Tuple(self + (dflt,), name=self.name)

    def _to7Tuple(self):
        '''(INTERNAL) see L{Rhumb8Tuple._to7Tuple}.
        '''
        return self


class Transform7Tuple(_Deprecated_NamedTuple):  # PYCHOK no cover
    '''DEPRECATED, use class L{Helmert7Tuple}, I{without} keyword arguments.'''
    _Names_ = _Helmert7Tuple._Names_
    _Units_ = _Helmert7Tuple._Units_

    def __new__(cls, tx=0, ty=0, tz=0, s=0,
                     sx=0, sy=0, sz=0, name=NN):
        t = map(_float, (tx, ty, tz, s, sx, sy, sz))
        return _Deprecated_NamedTuple.__new__(cls, *t, name=name)


class TriAngle4Tuple(_Deprecated_NamedTuple):
    '''DEPRECATED, use class L{TriAngle5Tuple}, ignoring item C{area}.'''
    assert _TriAngle5Tuple._Names_.index(_area_) == 4
    _Names_ = _TriAngle5Tuple._Names_[:4]
    _Units_ = _TriAngle5Tuple._Units_[:4]


class UtmUps4Tuple(_Deprecated_NamedTuple):  # PYCHOK no cover
    '''DEPRECATED and OBSOLETE, expect a L{UtmUps5Tuple} from method C{pygeodesy.Mgrs.toUtm(utm=None)}.

       4-Tuple C{(zone, hemipole, easting, northing)} with as C{zone} B{C{str}} and no C{band}.
    '''
    assert _UtmUps5Tuple._Names_.index(_band_) == 4
    _Names_ =                      _UtmUps5Tuple._Names_[ :4]  # band
    _Units_ = (_Deprecated_Str,) + _UtmUps5Tuple._Units_[1:4]


__all__ += _ALL_OTHER(ClipCS3Tuple, EasNorExact4Tuple, EcefCartesian,
                      HeightIDW, HeightIDW2, HeightIDW3, LatLonExact4Tuple,
                      Ned3Tuple, RefFrameError, Rhumb7Tuple, Transform7Tuple,
                      TriAngle4Tuple, UtmUps4Tuple)

# **) MIT License
#
# Copyright (C) 2018-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
