
# -*- coding: utf-8 -*-

u'''DEPRECATED constants, classes, functions, methods, etc.

Kept and exported for backward compatibility, including deprecated modules
C{pygeodesy.bases}, C{pygeodesy.datum} and C{pygeodesy,nvector}, previously
inside the C{pygeodesy} package.

Use either C{from pygeodesy import bases} or C{from pygeodesy.deprecated import
bases}.  Likewise for C{datum} and C{nvector}.
'''
from pygeodesy.interns import EPS, EPS_2, NAN, NN, R_M, _azi12_, _COMMASPACE_, \
                             _convergence_, _down_, _east_, _easting_, _end_, \
                             _float, _hemipole_, _lat_, _lat1_, _lat2_, _lon_, \
                             _lon1_, _lon2_, _negative_, _north_, _northing_, \
                             _s_, _s12_, _S12_, _scale_, _scalar_, _sep_, _SPACE_, \
                             _start_, _sx_, _sy_, _sz_, _tx_,  _ty_,  _tz_, \
                             _UNDER_, _value_, _zone_, _1_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, isLazy
from pygeodesy.named import _NamedTuple, _Pass
from pygeodesy.props import deprecated_class, deprecated_function, deprecated_method
from pygeodesy.units import Degrees, Easting, Float, Lat, Lon, Meter, Northing, \
                            Number_, Scalar, Scalar_, Str
if isLazy:  # XXX force import of all deprecated modules
    import pygeodesy.deprecated.bases as bases, \
           pygeodesy.deprecated.datum as datum, \
           pygeodesy.deprecated.nvector as nvector  # PYCHOK unused
    # XXX instead, use module_property or enhance .lazily

__all__ = _ALL_LAZY.deprecated
__version__ = '22.08.01'

EPS1_2 = _1_0 - EPS_2  # DEPRECATED
OK     = 'OK'          # DEPRECATED
_WGS84 = _UTM = object()


class _DeprecatedNamedTuple(_NamedTuple):
    '''(INTERNAL) Base class for DEPRECATED C{_NamedTuple} classes.
    '''
    def __new__(cls, *args, **kwds):
        deprecated_class(cls)
        return _NamedTuple.__new__(cls, *args, **kwds)


# DEPRECATED classes, for export and backward compatibility only
class ClipCS3Tuple(_DeprecatedNamedTuple):  # PYCHOK no cover
    '''DEPRECATED, see function L{pygeodesy.clipCS3}.'''
    _Names_ = (_start_, _end_, 'index')
    _Units_ = (_Pass,   _Pass,  Number_)


def EcefCartesian(*args, **kwds):
    '''DEPRECATED, use class L{LocalCartesian}.'''
    LocalCartesian = _MODS.ltp.LocalCartesian

    class EcefCartesian(LocalCartesian):
        '''DEPRECATED, use class L{LocalCartesian}.

           @note: This class is named I{incorrectly}, since it provides conversion to
                  and from I{local} cartesian coordinates in a I{local tangent plane}
                  and I{not geocentric} (ECEF) ones, as the name suggests.
        '''
        def __init__(self, latlonh0=0, lon0=0, height0=0, ecef=None, name=NN):
            deprecated_class(self.__class__)
            LocalCartesian.__init__(self, latlonh0=latlonh0, lon0=lon0, height0=height0, ecef=ecef, name=name)

        @deprecated_method
        def forward(self, latlonh, lon=None, height=0, M=False, name=NN):
            '''DEPRECATED, use method L{LocalCartesian.forward}.

               @return: I{Incorrectly}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C,
                        M, datum)} with I{local} C{(x, y, z)} coordinates for the given
                        I{geodetic} ones C{(lat, lon, height)}, case C{C=0} always,
                        optionally I{concatenated} L{EcefMatrix} C{M} and C{datum}.
            '''
            t = LocalCartesian.forward(self, latlonh, lon=lon, height=height, M=M, name=name)
            return _MODS.ecef.Ecef9Tuple(t.x, t.y, t.z, t.lat, t.lon, t.height,
                                                        0, t.M, t.ecef.datum,
                                                        name=t.name or self.name)

        @deprecated_method
        def reverse(self, xyz, y=None, z=None, M=False, name=NN):
            '''DEPRECATED, use method L{LocalCartesian.reverse}.

               @return: I{Incorrectly}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C,
                        M, datum)} with I{geodetic} coordinates C{(lat, lon, height)} for
                        the given I{local} ones C{(x, y, z)}, case C{C}, optionally
                        I{concatenated} L{EcefMatrix} C{M} and C{datum}.
            '''
            t = LocalCartesian.reverse(self, xyz, y=y, z=z, M=M, name=name)
            return _MODS.ecef.Ecef9Tuple(t.x, t.y, t.z, t.lat, t.lon, t.height,
                                                        t.ecef.C, t.M, t.ecef.datum,
                                                        name=t.name or self.name)

    _MODS.deprecated.EcefCartesian = EcefCartesian
    return EcefCartesian(*args, **kwds)


class EasNorExact4Tuple(_DeprecatedNamedTuple):
    '''DEPRECATED, use class L{Forward4Tuple}.'''
    _Names_ = (_easting_, _northing_, _convergence_, _scale_)
    _Units_ = ( Easting,   Northing,   Degrees,       Scalar)


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


class LatLonExact4Tuple(_DeprecatedNamedTuple):
    '''DEPRECATED, used class L{Reverse4Tuple}.'''
    _Names_ = (_lat_, _lon_, _convergence_, _scale_)
    _Units_ = ( Lat,   Lon,   Degrees,       Scalar)


class Ned3Tuple(_DeprecatedNamedTuple):  # was in .ellipsoidalNvector
    '''DEPRECATED, use class L{pygeodesy.Ned4Tuple}.'''
    _Names_ = (_north_, _east_,  _down_)
    _Units_ = ( Meter,   Meter,   Meter)


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


class Rhumb7Tuple(_DeprecatedNamedTuple):
    '''DEPRECATED, use class L{Rhumb8Tuple} ignoring item C{a12}.'''
    _Names_ = (_lat1_, _lon1_, _lat2_, _lon2_, _azi12_, _s12_, _S12_)
    _Units_ = (_Pass,  _Pass,  _Pass,  _Pass,  _Pass,   _Pass, _Pass)

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
        return _MODS.rhumbx.Rhumb8Tuple(self + (dflt,), name=self.name)

    def _to7Tuple(self, **unused):
        '''(INTERNAL) see L{Rhumb8Tuple._to7Tuple}.
        '''
        return self


class Transform7Tuple(_DeprecatedNamedTuple):  # PYCHOK no cover
    '''DEPRECATED, use class L{Helmert7Tuple} without keyword arguments.'''
    _Names_ = (_tx_,  _ty_,  _tz_,  _s_,   _sx_,  _sy_,  _sz_)
    _Units_ = ( Float, Float, Float, Float, Float, Float, Float)

    def __new__(cls, tx=0, ty=0, tz=0, s=0,
                     sx=0, sy=0, sz=0, name=NN):
        t = map(_float, (tx, ty, tz, s, sx, sy, sz))
        return _DeprecatedNamedTuple.__new__(cls, *t, name=name)


class UtmUps4Tuple(_DeprecatedNamedTuple):  # PYCHOK no cover
    '''OBSOLETE, expect a L{UtmUps5Tuple} from method C{Mgrs.toUtm(utm=None)}.

       4-Tuple C{(zone, hemipole, easting, northing)} as C{str},
       C{str}, C{meter} and C{meter}.
    '''
    _Names_ = (_zone_, _hemipole_, _easting_, _northing_)  # band
    _Units_ = ( Str,    Str,        Easting,   Northing)


@deprecated_function
def anStr(name, OKd='._-', sub=_UNDER_):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.anstr}.'''
    return _MODS.streprs.anstr(name, OKd=OKd, sub=sub)


@deprecated_function
def areaof(points, adjust=True, radius=R_M, wrap=True):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.areaOf}.'''
    return _MODS.points.areaOf(points, adjust=adjust, radius=radius, wrap=wrap)


@deprecated_function
def bounds(points, wrap=True, LatLon=None):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.boundsOf}.

       @return: 2-Tuple C{(latlonSW, latlonNE)} as B{C{LatLon}}
                or 4-Tuple C{(latS, lonW, latN, lonE)} if
                B{C{LatLon}} is C{None}.
    '''
    return tuple(_MODS.points.boundsOf(points, wrap=wrap, LatLon=LatLon))


@deprecated_function
def clipCS3(points, lowerleft, upperright, closed=False, inull=False):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.clipCS4}.

       @return: Yield a L{ClipCS3Tuple}C{(start, end, index)} for each
                edge of the I{clipped} path.
    '''
    for p1, p2, _, j in _MODS.clipy.clipCS4(points, lowerleft, upperright,
                                                    closed=closed, inull=inull):
        yield ClipCS3Tuple(p1, p2, j)


@deprecated_function
def clipDMS(deg, limit):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.clipDegrees}.'''
    return _MODS.dms.clipDegrees(deg, limit)


@deprecated_function
def clipStr(bstr, limit=50, white=NN):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.clips}.'''
    return _MODS.basics.clips(bstr, limit=limit, white=white)


@deprecated_function
def collins(pointA, pointB, pointC, alpha, beta, **useZ_Clas_and_kwds):
    '''DEPRECATED, use function L{pygeodesy.collins5}.'''
    from pygeodesy.resections import collins5  # _MODS.resections.collins5
    return collins5(pointA, pointB, pointC, alpha, beta, **useZ_Clas_and_kwds)


@deprecated_function
def copysign(x, y):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.copysign0}.'''
    return _MODS.basics.copysign0(x, y)


@deprecated_function
def decodeEPSG2(arg):  # PYCHOK no cover
    '''DEPRECATED, use function L{epsg.decode2}.

       @return: 2-Tuple C{(zone, hemipole)}
    '''
    return tuple(_MODS.epsg.decode2(arg))


@deprecated_function
def encodeEPSG(zone, hemipole=NN, band=NN):  # PYCHOK no cover
    '''DEPRECATED, use function L{epsg.encode}.

       @return: C{EPSG} code (C{int}).
    '''
    return int(_MODS.epsg.encode(zone, hemipole=hemipole, band=band))


@deprecated_function
def enStr2(easting, northing, prec, *extras):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.enstr2}.'''
    return _MODS.streprs.enstr2(easting, northing, (int(prec) // 2 - 5), *extras)


@deprecated_function
def equirectangular3(lat1, lon1, lat2, lon2, **options):  # PYCHOK no cover
    '''DEPRECATED, use function C{equirectangular_}.

       @return: 3-Tuple C{(distance2, delta_lat, delta_lon)}.
    '''
    return tuple(_MODS.formy.equirectangular_(lat1, lon1, lat2, lon2, **options)[:3])


@deprecated_function
def false2f(value, name=_value_, false=True, Error=ValueError):  # PYCHOK no cover
    '''DEPRECATED, use function L{falsed2f}.'''
    return falsed2f(falsed=false, Error=Error, **{name: value})


@deprecated_function
def falsed2f(falsed=True, Error=ValueError, **name_value):  # PYCHOK no cover
    '''DEPRECATED, use class L{Easting} or L{Northing}.

       Convert a falsed east-/northing to non-negative C{float}.

       @kwarg falsed: Value includes false origin (C{bool}).
       @kwarg Error: Optional, overriding error (C{Exception}).
       @kwarg name_value: One C{B{name}=value} pair.

       @return: The value (C{float}).

       @raise Error: Invalid or negative C{B{name}=value}.
    '''
    t = NN
    if len(name_value) == 1:
        try:
            for f in name_value.values():
                f = float(f)
                if falsed and f < 0:
                    break
                return f
            t = _COMMASPACE_('falsed', _negative_)
        except (TypeError, ValueError) as x:
            t = str(x)
    raise _MODS.errors._InvalidError(Error=Error, txt=t, **name_value)


@deprecated_function
def fStr(floats, prec=6, fmt=_MODS.streprs.Fmt.f, ints=False, sep=_COMMASPACE_):  # PYCHOK no cover
    '''DEPRECATED, use function L{fstr}.'''
    return _MODS.streprs.fstr(floats, prec=prec, fmt=fmt, ints=ints, sep=sep)


@deprecated_function
def fStrzs(floatstr):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.fstrzs}.'''
    return _MODS.streprs.fstrzs(floatstr)


@deprecated_function
def hypot3(x, y, z):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.hypot_}.'''
    return _MODS.fmath.hypot_(x, y, z)


@deprecated_function
def inStr(inst, *args, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.instr}.'''
    return _MODS.streprs.instr(inst, *args, **kwds)


@deprecated_function
def isenclosedby(point, points, wrap=False):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.isenclosedBy}.'''
    return _MODS.points.isenclosedBy(point, points, wrap=wrap)


@deprecated_function
def joined(*words, **sep):  # PYCHOK no cover
    '''DEPRECATED, use C{NN(...)}, C{NN.join_} or C{B{sep}.join}.'''
    return sep.get(_sep_, NN).join(map(str, words))


@deprecated_function
def joined_(*words, **sep):  # PYCHOK no cover
    '''DEPRECATED, use C{_SPACE_(...)}, C{_SPACE_.join_} or C{B{sep}.join}, sep=" ".'''
    return sep.get(_sep_, _SPACE_).join(map(str, words))


@deprecated_function
def nearestOn3(point, points, closed=False, wrap=False, **options):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.nearestOn5}.

       @return: 3-Tuple C{(lat, lon, distance)}
    '''
    return tuple(_MODS.points.nearestOn5(point, points, closed=closed, wrap=wrap, **options)[:3])


@deprecated_function
def nearestOn4(point, points, closed=False, wrap=False, **options):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.nearestOn5}.

       @return: 4-Tuple C{(lat, lon, distance, angle)}
    '''
    return tuple(_MODS.points.nearestOn5(point, points, closed=closed, wrap=wrap, **options)[:4])


@deprecated_function
def parseUTM(strUTM, datum=_WGS84, Utm=_UTM, name=NN):  # PYCHOK no cover
    '''DEPRECATED, use function L{parseUTM5}.

       @return: The UTM coordinate (B{L{Utm}}) or 4-tuple C{(zone,
                hemisphere, easting, northing)} if B{C{Utm}} is C{None}.
    '''
    d = _MODS.datums.Datums.WGS84 if datum is _WGS84 else datum  # PYCHOK shadows?
    U = _MODS.utm.Utm if Utm is _UTM else Utm
    r = _MODS.utm.parseUTM5(strUTM, datum=d, Utm=U, name=name)
    if isinstance(r, tuple):  # UtmUps5Tuple
        r = r.zone, r.hemipole, r.easting, r.northing  # no band
    return r


@deprecated_function
def perimeterof(points, closed=False, adjust=True, radius=R_M, wrap=True):  # PYCHOK no cover
    '''DEPRECATED, use function L{perimeterOf}.'''
    return _MODS.points.perimeterOf(points, closed=closed, adjust=adjust, radius=radius, wrap=wrap)


@deprecated_function
def polygon(points, closed=True, base=None):  # PYCHOK no cover
    '''DEPRECATED, use function L{points2}.'''
    from pygeodesy.deprecated.bases import points2  # PYCHOK import
    return points2(points, closed=closed, base=base)


@deprecated_function
def scalar(value, low=EPS, high=1.0, name=_scalar_, Error=ValueError):  # PYCHOK no cover
    '''DEPRECATED, use class L{Number_} or L{Scalar_}.

       @return: New value (C{float} or C{int} for C{int} B{C{low}}).

       @raise Error: Invalid B{C{value}}.
    '''
    C_ = Number_ if _MODS.basics.isint(low) else Scalar_
    return C_(value, name=name, Error=Error, low=low, high=high)


@deprecated_function
def simplify2(points, pipe, radius=R_M, shortest=False, indices=False, **options):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.simplifyRW}.
    '''
    return _MODS.simplify.simplifyRW(points, pipe, radius=radius, shortest=shortest,
                                                   indices=indices, **options)


@deprecated_function
def tienstra(pointA, pointB, pointC, alpha, **beta_gamma_useZ_Clas_and_kwds):
    '''DEPRECATED, use function L{pygeodesy.tienstra7}.'''
    return _MODS.resections.tienstra7(pointA, pointB, pointC, alpha, **beta_gamma_useZ_Clas_and_kwds)


@deprecated_function
def toUtm(latlon, lon=None, datum=None, Utm=_UTM, cmoff=True, name=NN):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.toUtm8}.

       @return: The UTM coordinate (B{C{Utm}}) or a 6-tuple C{(zone,
                easting, northing, band, convergence, scale)} if
                B{C{Utm}} is C{None} or B{C{cmoff}} is C{False}.
    '''
    U = _MODS.utm.Utm if Utm is _UTM else Utm
    r = _MODS.utm.toUtm8(latlon, lon=lon, datum=datum, Utm=U, name=name, falsed=cmoff)
    if isinstance(r, tuple):  # UtmUps8Tuple
        # no hemisphere/pole and datum
        r = r.zone, r.easting, r.northing, r.band, r.convergence, r.scale
    return r


@deprecated_function
def unsign0(x):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.unsigned0}.'''
    return _MODS.basics.unsigned0(x)


@deprecated_function
def unStr(name, *args, **kwds):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.unstr}.'''
    return _MODS.streprs.unstr(name, *args, **kwds)


@deprecated_function
def utmZoneBand2(lat, lon):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.utmZoneBand5}.

       @return: 2-Tuple C{(zone, band)}.
    '''
    r = _MODS.utm.utmZoneBand5(lat, lon)  # UtmUpsLatLon5Tuple
    return r.zone, r.band

# **) MIT License
#
# Copyright (C) 2018-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
