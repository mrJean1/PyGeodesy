
# -*- coding: utf-8 -*-

u'''Named, I{Local Tangent Plane} (LTP) tuples.

Local coordinate classes L{XyzLocal}, L{Enu}, L{Ned} and L{Aer}
and local coordinate tuples L{Local9Tuple}, L{Xyz4Tuple}, L{Enu4Tuple},
L{Ned4Tuple}, L{Aer4Tuple} and L{Footprint5Tuple}.

@see: References in module L{ltp}.
'''

from pygeodesy.basics import issubclassof, _xinstanceof
from pygeodesy.dms import F_D, toDMS
from pygeodesy.errors import _TypesError, _xkwds
from pygeodesy.fmath import fdot, hypot, hypot_  # PYCHOK shared
from pygeodesy.interns import NN, _4_, _azimuth_, _COMMASPACE_, \
                             _down_, _east_, _ecef_, _elevation_, \
                             _height_, _lat_, _lon_, _ltp_, _M_, \
                             _name_, _north_, _up_, _x_, _xyz_, \
                             _y_, _z_, _0_0, _90_0, _N_90_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass, _xnamed
from pygeodesy.namedTuples import LatLon2Tuple, PhiLam2Tuple, Vector3Tuple
from pygeodesy.props import deprecated_class, deprecated_method, \
                            deprecated_Property_RO, Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, strs, _xzipairs
from pygeodesy.units import Bearing, Degrees, Degrees_, Height, Lat, \
                            Lon, Meter, Meter_
from pygeodesy.utily import atan2d, atan2b, sincos2d_
from pygeodesy.vector3d import Vector3d

from math import cos, radians

__all__ = _ALL_LAZY.ltpTuples
__version__ = '21.12.28'

_aer_        = 'aer'
_alt_        = 'alt'
_roll_       = 'roll'
_slantrange_ = 'slantrange'
_tilt_       = 'tilt'
_yaw_        = 'yaw'


def _er2gr(e, r):
    '''(INTERNAL) Elevation and slant range to ground range.
    '''
    c = cos(radians(e))
    return Meter_(groundrange=r * c)


def _4Tuple2Cls(inst, Cls, Cls_kwds):
    '''(INTERNAL) Convert 4-Tuple to C{Cls} instance.
    '''
    if Cls is None:
        return inst
    elif issubclassof(Cls, Aer):
        return inst.xyzLocal.toAer(Aer=Cls, **Cls_kwds)
    elif issubclassof(Cls, Enu):
        return inst.xyzLocal.toEnu(Enu=Cls, **Cls_kwds)
    elif issubclassof(Cls, Ned):
        return inst.xyzLocal.toNed(Ned=Cls, **Cls_kwds)
    elif issubclassof(Cls, XyzLocal):
        return inst.xyzLocal.toXyz(Xyz=Cls, **Cls_kwds)
    elif Cls is Local9Tuple:
        return inst.xyzLocal.toLocal9Tuple(**Cls_kwds)
    n = inst.__class__.__name__[:3]
    raise _TypesError(n, Cls, Aer, Enu, Ned, XyzLocal)


def _xyz2aer4(inst):
    '''(INTERNAL) Convert C{(x, y, z}) to C{(A, E, R)}.
    '''
    x, y, z, _ = inst.xyz4
    A = Bearing(azimuth=atan2b(x, y))
    E = Degrees(elevation=atan2d(z, hypot(x, y)))
    R = Meter(slantrange=hypot_(x, y, z))
    return Aer4Tuple(A, E, R, inst.ltp, name=inst.name)


class Aer(_NamedBase):
    '''Local C{Azimuth-Elevation-Range} (AER) in a I{local tangent plane}.
    '''
    _azimuth     = _0_0   # bearing from North (C{degrees360})
    _elevation   = _0_0   # tilt, pitch from horizon (C{degrees}).
    _ltp         =  None  # local tangent plane (C{Ltp}), origin
    _slantrange  = _0_0   # distance (C{Meter})
    _toStr       = _aer_

    def __init__(self, aer, elevation=_0_0, slantrange=_0_0, ltp=None, name=NN):
        '''New L{Aer}.

           @arg aer: Scalar azimuth, bearing from North (C{degrees360}) or
                     local (L{Aer}, L{Aer4Tuple}).
           @kwarg elevation: Scalar angle I{above} horizon, I{above} B{C{ltp}}
                             (C{degrees90}) if scalar B{C{aer}}.
           @kwarg slantrange: Scalar distance (C{meter}) if scalar B{C{aer}}.
           @kwarg ltp: The I{local tangent plane}, (geodetic) origin (L{Ltp},
                       L{LocalCartesian}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{ltp}}.

           @raise UnitError: Invalid B{C{azimuth}}, B{C{elevation}} or
                             or B{C{slantrange}}.
        '''
        try:
            self._azimuth, self._elevation, self._slantrange = \
              aer.azimuth,   aer.elevation,   aer.slantrange
            _xinstanceof(Aer, Aer4Tuple, aer=aer)
            p = getattr(aer, _ltp_, ltp)
            n = getattr(aer, _name_, name)
        except AttributeError:
            self._azimuth    = Bearing(azimuth=aer)
            self._elevation  = Degrees_(elevation=elevation, low=_N_90_0, high=_90_0)
            self._slantrange = Meter_(slantrange=slantrange)
            p, n = ltp, name

        if p:
            self._ltp = _MODS.ltp._xLtp(p)
        if name:
            self.name = n

    @Property_RO
    def aer4(self):
        '''Get the C{(azimuth, elevation, slantrange, ltp)} components (L{Aer4Tuple}).
        '''
        return Aer4Tuple(self._azimuth, self._elevation, self._slantrange, self.ltp, name=self.name)

    @Property_RO
    def azimuth(self):
        '''Get the Azimuth, bearing from North (C{degrees360}).
        '''
        return self._azimuth

    @Property_RO
    def down(self):
        '''Get the Down component (C{meter}).
        '''
        return self.xyzLocal.down

    @Property_RO
    def east(self):
        '''Get the East component (C{meter}).
        '''
        return self.xyzLocal.east

    @Property_RO
    def elevation(self):
        '''Get the Elevation, tilt above horizon (C{degrees90}).
        '''
        return self._elevation

    @Property_RO
    def groundrange(self):
        '''Get the I{ground range}, distance (C{meter}).
        '''
        return _er2gr(self._elevation, self._slantrange)

    @Property_RO
    def ltp(self):
        '''Get the I{local tangent plane} (L{Ltp}).
        '''
        return self._ltp

    @Property_RO
    def north(self):
        '''Get the North component (C{meter}).
        '''
        return self.xyzLocal.north

    @Property_RO
    def slantrange(self):
        '''Get the I{slant Range}, distance (C{meter}).
        '''
        return self._slantrange

    def toRepr(self, prec=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this AER as azimuth
           (bearing), elevation and slant range.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between AERs (C{str}).

           @return: This AER as "[A:degrees360, E:degrees90, R:meter]" (C{str}).
        '''
        t = (toDMS(self.azimuth,    form=F_D, prec=prec, ddd=0),
             toDMS(self.elevation,  form=F_D, prec=prec, ddd=0),
             fstr( self.slantrange, prec=3 if prec is None else prec))
        return _xzipairs(self._toStr.upper(), t, sep=sep, fmt=fmt)

    def toStr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this ENU/NED/AER.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between NEDs (C{str}).

           @return: This AER as "[degrees360, degrees90, meter]" (C{str}).
        '''
        a = self._toStr
        t = getattr(self, a + _4_)[:len(a)]
        return fmt(sep.join(strs(t, prec=prec)))

    def toXyz(self, Xyz=None, **Xyz_kwds):
        '''Get the local I{X, Y, Z} (XYZ) components.

           @kwarg Xyz: Class to return XYZ (L{XyzLocal}, L{Enu},
                       L{Ned}, L{Aer}) or C{None}.
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: XYZ as an B{C{Xyz}} instance or if C{B{Xyz} is None},
                    an L{Xyz4Tuple}C{(x, y, z, ltp)}.

           @raise TypeError: Invalid B{C{Xyz}}.
        '''
        return self.xyz4._toXyz(Xyz, Xyz_kwds)

    @Property_RO
    def up(self):
        '''Get the Up component (C{meter}).
        '''
        return self.xyzLocal.up

    @Property_RO
    def x(self):
        '''Get the X component (C{meter}).
        '''
        return self.xyz4.x

    @Property_RO
    def xyz(self):
        '''Get the I{local} C{(X, Y, Z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)  # like Ecef9Tuple.xyz, Local6tuple.xyz

    @Property_RO
    def xyz4(self):
        '''Get the C{(x, y, z, ltp)} components (L{Xyz4Tuple}).
        '''
        sA, cA, sE, cE = sincos2d_(self._azimuth, self._elevation)
        R = self._slantrange
        r = cE * R  # ground range
        return Xyz4Tuple(sA * r, cA * r, sE * R, self.ltp, name=self.name)

    @Property_RO
    def xyzLocal(self):
        '''Get this AER as an L{XyzLocal}.
        '''
        return XyzLocal(self.xyz4)

    @Property_RO
    def y(self):
        '''Get the Y component (C{meter}).
        '''
        return self.xyz4.y

    @Property_RO
    def z(self):
        '''Get the Z component (C{meter}).
        '''
        return self.xyz4.z


class Aer4Tuple(_NamedTuple):
    '''4-Tuple C{(azimuth, elevation, slantrange, ltp)},
       all in C{meter} except C{ltp}.
    '''
    _Names_ = (_azimuth_, _elevation_, _slantrange_, _ltp_)
    _Units_ = ( Meter,     Meter,       Meter,       _Pass)

    def _toAer(self, Cls, Cls_kwds):
        '''(INTERNAL) Return C{Cls(..., **Cls_kwds)} instance.
        '''
        if issubclassof(Cls, Aer):
            return Cls(*self, **_xkwds(Cls_kwds, name=self.name))
        else:
            return _4Tuple2Cls(self, Cls, Cls_kwds)

    @Property_RO
    def groundrange(self):
        '''Get the I{ground range}, distance (C{meter}).
        '''
        return _er2gr(self.elevation, self.slantrange)  # PYCHOK _Tuple

    @Property_RO
    def xyzLocal(self):
        '''Get this L{Aer4Tuple} as an L{XyzLocal}.
        '''
        return Aer(self).xyzLocal


class Attitude4Tuple(_NamedTuple):
    '''4-Tuple C{(alt, tilt, yaw, roll)} with C{altitude} in (positive)
       C{meter} and C{tilt}, C{yaw} and C{roll} in C{degrees} representing
       the attitude of a plane or camera.
    '''
    _Names_ = (_alt_, _tilt_,  _yaw_,   _roll_)
    _Units_ = ( Meter, Bearing, Degrees, Degrees)

    @Property_RO
    def atyr(self):
        '''Return this attitude (L{Attitude4Tuple}).
        '''
        return self

    @Property_RO
    def tyr3d(self):
        '''Get this attitude's (3-D) directional vector (L{Vector3d}).
        '''
        return _MODS.ltp.Attitude(self).tyr3d


class Ned(_NamedBase):
    '''Local C{North-Eeast-Down} (NED) location in a I{local tangent plane}.

       @see: L{Enu} and L{Ltp}.
    '''
    _down  = _0_0   # down, -XyzLocal.z (C{meter}).
    _east  = _0_0   # east, XyzLocal.y (C{meter}).
    _ltp   =  None  # local tangent plane (C{Ltp}), origin
    _north = _0_0   # north, XyzLocal.x (C{meter})
    _toStr = 'ned'

    def __init__(self, ned, east=_0_0, down=_0_0, ltp=None, name=NN):
        '''New L{Ned} vector.

           @arg ned: Scalar north component (C{meter}) or local NED (L{Ned},
                     L{Ned4Tuple}).
           @kwarg east: Scalar east component (C{meter}) if scalar B{C{ned}}.
           @kwarg down: Scalar down component, normal to I{inside} surface of
                        the ellipsoid or sphere (C{meter}) if scalar B{C{ned}}.
           @kwarg ltp: The I{local tangent plane}, (geodetic) origin (L{Ltp},
                       L{LocalCartesian}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{ltp}}.

           @raise UnitError: Invalid B{C{north}}, B{C{east}} or B{C{down}}.
        '''
        try:
            self._north, self._east, self._down = ned.north, ned.east, ned.down
            _xinstanceof(Ned, Ned4Tuple, ned=ned)
            p = getattr(ned, _ltp_, ltp)
            n = getattr(ned, _name_, name)
        except AttributeError:
            self._north = Meter(north=ned or _0_0)
            self._east  = Meter(east=east or _0_0)
            self._down  = Meter(down=down or _0_0)
            p, n = ltp, name

        if p:
            self._ltp = _MODS.ltp._xLtp(p)
        if n:
            self.name = n

    @Property_RO
    def aer4(self):
        '''Get the C{(azimuth, elevation, slantrange, ltp)} components (L{Aer4Tuple}).
        '''
        return _xyz2aer4(self)

    @Property_RO
    def azimuth(self):
        '''Get the Azimuth, bearing from North (C{degrees360}).
        '''
        return self.aer4.azimuth

    @deprecated_Property_RO
    def bearing(self):
        '''DEPRECATED, use C{azimuth}.'''
        return self.azimuth

    @Property_RO
    def down(self):
        '''Get the Down component (C{meter}).
        '''
        return self._down

    @Property_RO
    def east(self):
        '''Get the East component (C{meter}).
        '''
        return self._east

    @Property_RO
    def elevation(self):
        '''Get the Elevation, tilt above horizon (C{degrees90}).
        '''
        return self.aer4.elevation  # neg(degrees90(asin1(self.down / self.length))))

    @Property_RO
    def groundrange(self):
        '''Get the I{ground range}, distance (C{meter}).
        '''
        return Meter(groundrange=hypot(self.north, self.east))

    @deprecated_Property_RO
    def length(self):
        '''DEPRECATED, use C{slantrange}.'''
        return self.slantrange

    @Property_RO
    def ltp(self):
        '''Get the I{local tangent plane} (L{Ltp}).
        '''
        return self._ltp

    @deprecated_Property_RO
    def ned(self):
        '''DEPRECATED, use property C{ned4}.'''
        return Ned3Tuple(self.north, self.east, self.down, name=self.name)

    @Property_RO
    def ned4(self):
        '''Get the C{(north, east, down, ltp)} components (L{Ned4Tuple}).
        '''
        return Ned4Tuple(self.north, self.east, self.down, self.ltp, name=self.name)

    @Property_RO
    def north(self):
        '''Get the North component (C{meter}).
        '''
        return self._north

    @Property_RO
    def slantrange(self):
        '''Get the I{slant Range}, distance (C{meter}).
        '''
        return self.aer4.slantrange

    @deprecated_method
    def to3ned(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{ned4}.'''
        return self.ned  # XXX deprecated too

    def toRepr(self, prec=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this NED.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between items (C{str}).

           @return: This NED as "[N:meter, E:meter, D:meter]" (C{str}).
        '''
        a = self._toStr
        t = getattr(self, a + _4_)[:len(a)]
        t = strs(t, prec=3 if prec is None else prec)
        return _xzipairs(a.upper(), t, sep=sep, fmt=fmt)

    def toStr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this NED.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between NEDs (C{str}).

           @return: This NED as "[meter, meter, meter]" (C{str}).
        '''
        a = self._toStr
        t = getattr(self, a + _4_)[:len(a)]
        return fmt(sep.join(strs(t, prec=prec)))

    @deprecated_method
    def toVector3d(self):
        '''DEPRECATED, use property L{xyz}.'''
        return self.xyz

    def toXyz(self, Xyz=None, **Xyz_kwds):
        '''Get the local I{X, Y, Z} (XYZ) components.

           @kwarg Xyz: Class to return XYZ (L{XyzLocal}, L{Enu},
                       L{Ned}, L{Aer}) or C{None}.
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: XYZ as an B{C{Xyz}} instance or if C{B{Xyz} is None},
                    an L{Xyz4Tuple}C{(x, y, z, ltp)}.

           @raise TypeError: Invalid B{C{Xyz}}.
        '''
        return self.xyz4._toXyz(Xyz, Xyz_kwds)

    @Property_RO
    def up(self):
        '''Get the Up component (C{meter}).
        '''
        return Meter(up=-self._down)  # negated

    @Property_RO
    def x(self):
        '''Get the X component (C{meter}).
        '''
        return Meter(x=self._east)  # 2nd arg, E

    @Property_RO
    def xyz(self):
        '''Get the I{local} C{(X, Y, Z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)  # like Ecef9Tuple.xyz, Local6tuple.xyz

    @Property_RO
    def xyz4(self):
        '''Get the C{(x, y, z, ltp)} components (L{Xyz4Tuple}).
        '''
        return Xyz4Tuple(self.x, self.y, self.z, self.ltp, name=self.name)

    @Property_RO
    def xyzLocal(self):
        '''Get this NED as an L{XyzLocal}.
        '''
        return XyzLocal(self.xyz4)

    @Property_RO
    def y(self):
        '''Get the Y component (C{meter}).
        '''
        return Meter(y=self._north)  # 1st arg N

    @Property_RO
    def z(self):
        '''Get the Z component (C{meter}).
        '''
        return Meter(z=-self._down)  # negated


class Ned3Tuple(_NamedTuple):  # was in .ellipsoidalNvector
    '''3-Tuple C{(north, east, down)}.  DEPRECATED, use L{pygeodesy.Ned4Tuple}.
    '''
    _Names_ = (_north_, _east_,  _down_)
    _Units_ = ( Meter,   Meter,   Meter)

    def __new__(cls, north, east, down, name=NN):
        deprecated_class(cls)
        return _NamedTuple.__new__(cls, north, east, down, name=name)


class Ned4Tuple(_NamedTuple):
    '''4-Tuple C{(north, east, down, ltp)}, all in C{meter} except C{ltp}.
    '''
    _Names_ = (_north_, _east_, _down_, _ltp_)
    _Units_ = ( Meter,   Meter,  Meter, _Pass)

    def _toNed(self, Cls, Cls_kwds):
        '''(INTERNAL) Return C{Cls(..., **Cls_kwds)} instance.
        '''
        if issubclassof(Cls, Ned):
            return Cls(*self, **_xkwds(Cls_kwds, name=self.name))
        else:
            return _4Tuple2Cls(self, Cls, Cls_kwds)

    @Property_RO
    def xyzLocal(self):
        '''Get this L{Ned4Tuple} as an L{XyzLocal}.
        '''
        return Ned(self).xyzLocal


class XyzLocal(Vector3d):
    '''Local C{(x, y, z)} in a I{local tangent plane} (LTP)
       and base class for local L{Enu}, L{Ned} and L{Aer}.
    '''
    _ltp   =  None  # local tangent plane (C{Ltp}), origin
    _toStr = _xyz_

    def __init__(self, xyz, y=_0_0, z=_0_0, ltp=None, name=NN):
        '''New L{XyzLocal}.

           @arg xyz: Scalar X component (C{meter}), C{east} or a local (L{XyzLocal},
                     L{Xyz4Tuple}, L{Enu}, L{Enu4Tuple}, L{Local9Tuple}).
           @kwarg y: Scalar Y component (C{meter}) if scalar B{C{xyz}}, C{north}.
           @kwarg z: Scalar Z component if scalar B{C{xyz}}, normal C{up} from the
                     surface of the ellipsoid or sphere (C{meter}).
           @kwarg ltp: The I{local tangent plane}, (geodetic) origin (L{Ltp},
                       L{LocalCartesian}).

           @raise TypeError: Invalid B{C{xyz}}, B{C{ltp}}.

           @raise UnitError: Invalid B{C{x}}, B{C{y}} or B{C{z}}.
        '''
        try:
            self._x, self._y, self._z = xyz.x, xyz.y, xyz.z
            _xinstanceof(XyzLocal, Xyz4Tuple, Enu, Enu4Tuple, Local9Tuple, Ned,
                      **{self._toStr: xyz})
            p = getattr(xyz, _ltp_, ltp)
            n = getattr(xyz, _name_, name)
        except AttributeError:
            self._x = Meter(x=xyz or _0_0)
            self._y = Meter(y=y or _0_0)
            self._z = Meter(z=z or _0_0)
            p, n = ltp, name

        if p:
            self._ltp = _MODS.ltp._xLtp(p)
        if n:
            self.name = n

    def __str__(self):
        return self.toStr()

    @Property_RO
    def aer4(self):
        '''Get the C{(azimuth, elevation, slantrange, ltp)} components (L{Aer4Tuple}).
        '''
        return _xyz2aer4(self)

    @Property_RO
    def azimuth(self):
        '''Get the Azimuth, bearing from North (C{degrees360}).

           @see: U{Azimuth<https://GSSC.ESA.int/navipedia/index.php/
                 Transformations_between_ECEF_and_ENU_coordinates>}.
        '''
        return self.aer4.azimuth

    def classof(self, *args, **kwds):
        '''Create another instance of this very class.

           @arg args: Optional, positional arguments.
           @kwarg kwds: Optional, keyword arguments.

           @return: New instance (C{self.__class__}).
        '''
        kwds = _xkwds(kwds, ltp=self.ltp, name=self.name)
        return self.__class__(*args, **kwds)

    @Property_RO
    def down(self):
        '''Get the Down component (C{meter}).
        '''
        return Meter(down=-self.z)

    @property_RO
    def ecef(self):
        '''Get this LTP's ECEF converter (C{Ecef...} I{instance}).
        '''
        return self.ltp.ecef

    @Property_RO
    def east(self):
        '''Get the East component (C{meter}).
        '''
        return Meter(east=self.x)

    @Property_RO
    def elevation(self):
        '''Get the Elevation, tilt above horizon (C{degrees90}).

           @see: U{Elevation<https://GSSC.ESA.int/navipedia/index.php/
                 Transformations_between_ECEF_and_ENU_coordinates>}.
        '''
        return self.aer4.elevation  # neg(degrees90(asin1(self.down / self.length))))

    @Property_RO
    def enu4(self):
        '''Get the C{(east, north, up, ltp)} components (L{Enu4Tuple}).
        '''
        return Enu4Tuple(self.east, self.north, self.up, self.ltp, name=self.name)

    @Property_RO
    def groundrange(self):
        '''Get the I{ground range}, distance (C{meter}).
        '''
        return Meter(groundrange=hypot(self.x, self.y))

    @Property_RO
    def ltp(self):
        '''Get the I{local tangent plane} (L{Ltp}).
        '''
        return self._ltp

    @Property_RO
    def ned4(self):
        '''Get the C{(north, east, down, ltp)} components (L{Ned4Tuple}).
        '''
        return Ned4Tuple(self.north, self.east, self.down, self.ltp, name=self.name)

    @Property_RO
    def north(self):
        '''Get the North component (C{meter}).
        '''
        return Meter(north=self.y)

    @Property_RO
    def slantrange(self):
        '''Get the I{slant Range}, distance (C{meter}).
        '''
        return self.aer4.slantrange

    def toAer(self, Aer=None, **Aer_kwds):
        '''Get the local I{Azimuth, Elevation, slantRange} components.

           @kwarg Aer: Class to return AER (L{Aer}) or C{None}.
           @kwarg Aer_kwds: Optional, additional B{C{Aer}} keyword
                            arguments, ignored if C{B{Aer} is None}.

           @return: AER as an L{Aer} instance or if C{B{Aer} is None}, an
                    L{Aer4Tuple}C{(azimuth, elevation, slantrange, ltp)}.

           @raise TypeError: Invalid B{C{Aer}}.
        '''
        return self.aer4._toAer(Aer, Aer_kwds)

    def toCartesian(self, Cartesian=None, ltp=None, **Cartesian_kwds):
        '''Get the geocentric C{(x, y, z)} (ECEF) coordinates of this local.

           @kwarg Cartesian: Optional class to return C{(x, y, z)} (C{Cartesian})
                             or C{None}.
           @kwarg ltp: Optional I{local tangent plane} (LTP) (L{Ltp}),
                       overriding this C{ltp}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}} keyword
                                  arguments, ignored if C{B{Cartesian} is None}.

           @return: A B{C{Cartesian}} instance of if C{B{Cartesian} is None}, an
                    L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{M=None}, always.

           @raise TypeError: Invalid B{C{ltp}}, B{C{Cartesian}} or
                             B{C{Cartesian_kwds}} argument.
        '''
        ltp = _MODS.ltp._xLtp(self.ltp if ltp is None else ltp)
        if Cartesian is None:
            r = ltp._local2ecef(self, nine=True)
        else:
            x, y, z = ltp._local2ecef(self)
            kwds = _xkwds(Cartesian_kwds, datum=ltp.datum)
            r = Cartesian(x, y, z, **kwds)
        return _xnamed(r, self.name or ltp.name)

    def toEnu(self, Enu=None, **Enu_kwds):
        '''Get the local I{East, North, Up} (ENU) components.

           @kwarg Enu: Class to return ENU (L{Enu}) or C{None}.
           @kwarg Enu_kwds: Optional, additional B{C{Enu}} keyword
                            arguments, ignored if C{B{Enu} is None}.

           @return: ENU as an L{Enu} instance or if C{B{Enu} is None},
                    an L{Enu4Tuple}C{(east, north, up, ltp)}.
        '''
        return self.enu4._toEnu(Enu, Enu_kwds)

    def toLatLon(self, LatLon=None, ltp=None, **LatLon_kwds):
        '''Get the geodetic C{(lat, lon, height)} coordinates if this local.

           @kwarg LatLon: Optional class to return C{(x, y, z)} (C{LatLon})
                          or C{None}.
           @kwarg ltp: Optional I{local tangent plane} (LTP) (L{Ltp}),
                       overriding this ENU/NED/AER/XYZ's LTP.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: An B{C{LatLon}} instance of if C{B{LatLon} is None}, an
                    L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{M=None}, always.

           @raise TypeError: Invalid B{C{ltp}}, B{C{LatLon}} or
                             B{C{LatLon_kwds}} argument.
        '''
        ltp = _MODS.ltp._xLtp(self.ltp if ltp is None else ltp)
        r = ltp._local2ecef(self, nine=True)
        if LatLon is None:
            r = _xnamed(r, self.name or ltp.name)
        else:
            kwds = _xkwds(LatLon_kwds, height=r.height, datum=r.datum,
                                       name=self.name or ltp.name)
            r = LatLon(r.lat, r.lon, **kwds)  # XXX ltp?
        return r

    def toLocal9Tuple(self, M=False, name=NN):
        '''Get this local as a C{Local9Tuple}.

           @kwarg M: Optionally include the rotation matrix (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: L{Local9Tuple}C{(x, y, z, lat, lon, height, ltp,
                    ecef, M)} with C{ltp} this C{Ltp}, C{ecef} an
                    L{Ecef9Tuple} and C{M} L{EcefMatrix} or C{None}.
        '''
        ltp = self.ltp  # see C{self.toLatLon}
        t = ltp._local2ecef(self, nine=True, M=M)
        return Local9Tuple(self.x, self.y, self.z, t.lat, t.lon, t.height,
                                           ltp, t, t.M, name=name or t.name)

    def toNed(self, Ned=None, **Ned_kwds):
        '''Get the local I{North, East, Down} (Ned) components.

           @kwarg Ned: Class to return NED (L{Ned}) or C{None}.
           @kwarg Ned_kwds: Optional, additional B{C{Ned}} keyword
                            arguments, ignored if C{B{Ned} is None}.

           @return: NED as an L{Ned} instance or if C{B{Ned} is None},
                    an L{Ned4Tuple}C{(north, east, down, ltp)}.
        '''
        return self.ned4._toNed(Ned, Ned_kwds)

    def toRepr(self, prec=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this ENU/NED/XYZ.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between items (C{str}).

           @return: This XYZ/ENU as "[E:meter, N:meter, U:meter]",
                    "[N:meter, E:meter, D:meter]" respectively
                    "[X:meter, Y:meter, Z:meter]" (C{str}).
        '''
        a = self._toStr
        t = getattr(self, a + _4_)[:len(a)]
        t = strs(t, prec=3 if prec is None else prec)
        return _xzipairs(a.upper(), t, sep=sep, fmt=fmt)

    def toStr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this ENU/XYZ.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between items (C{str}).

           @return: This XYZ/ENU as "[meter, meter, meter]" (C{str}).
        '''
        a = self._toStr
        t = getattr(self, a + _4_)[:len(a)]
        return fmt(sep.join(strs(t, prec=prec)))

    def toXyz(self, Xyz=None, **Xyz_kwds):
        '''Get the local I{X, Y, Z} (XYZ) components.

           @kwarg Xyz: Class to return XYZ (L{XyzLocal}, L{Enu},
                       L{Ned}, L{Aer}) or C{None}.
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: XYZ as an B{C{Xyz}} instance or if C{B{Xyz} is None},
                    an L{Xyz4Tuple}C{(x, y, z, ltp)}.
        '''
        return self.xyz4._toXyz(Xyz, Xyz_kwds)

    @Property_RO
    def up(self):
        '''Get the Up component (C{meter}).
        '''
        return Meter(up=self.z)

#   @Property_RO
#   def x(self):  # see: Vector3d.x
#       '''Get the X component (C{meter}).
#       '''
#       return self._x

#   @Property_RO
#   def xyz(self):  # see: Vector3d.xyz
#       '''Get the I{local} C{(X, Y, Z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
#       '''
#       return Vector3Tuple(self.x, self.y, self.z, name=self.name)  # like Ecef9Tuple.xyz, Local6tuple.xyz

    @Property_RO
    def xyz4(self):
        '''Get the C{(x, y, z, ltp)} components (L{Xyz4Tuple}).
        '''
        return Xyz4Tuple(self.x, self.y, self.z, self.ltp, name=self.name)

    @Property_RO
    def xyzLocal(self):
        '''Get this L{XyzLocal}.
        '''
        return self

#   @Property_RO
#   def y(self):  # see: Vector3d.y
#       '''Get the Y component (C{meter}).
#       '''
#       return self._y
#
#   @Property_RO
#   def z(self):  # see: Vector3d.z
#       '''Get the Z component (C{meter}).
#       '''
#       return self._z


class Xyz4Tuple(_NamedTuple):
    '''4-Tuple C{(x, y, z, ltp)}, all in C{meter} except C{ltp}.
    '''
    _Names_ = (_x_,   _y_,   _z_,    _ltp_)
    _Units_ = ( Meter, Meter, Meter, _Pass)

    def _toXyz(self, Cls, Cls_kwds):
        '''(INTERNAL) Return C{Cls(..., **Cls_kwds)} instance.
        '''
        if issubclassof(Cls, XyzLocal):
            return Cls(*self, **_xkwds(Cls_kwds, name=self.name))
        else:
            return _4Tuple2Cls(self, Cls, Cls_kwds)

    @Property_RO
    def xyzLocal(self):
        '''Get this L{Xyz4Tuple} as an L{XyzLocal}.
        '''
        return XyzLocal(self)


class Enu(XyzLocal):
    '''Local C{Eeast-North-Up} (ENU) location in a I{local tangent plane}.

       @see: U{East, North, Up (ENU)<https://GSSC.ESA.int/navipedia/index.php/
             Transformations_between_ECEF_and_ENU_coordinates>} coordinates.
    '''
    _toStr = 'enu'

    def __init__(self, enu, north=_0_0, up=_0_0, ltp=None, name=NN):
        '''New L{Enu}.

           @arg enu: Scalar East component (C{meter}) or a local (L{Enu},
                     (L{XyzLocal}, L{Local9Tuple}).
           @kwarg north: Scalar North component (C{meter}) if scalar B{C{enu}}.
           @kwarg up: Scalar Up component if scalar B{C{enu}}, normal from the
                      surface of the ellipsoid or sphere (C{meter}).
           @kwarg ltp: The I{local tangent plane}, (geodetic) origin (L{Ltp},
                       L{LocalCartesian}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{enu}} or B{C{ltp}}.

           @raise UnitError: Invalid B{C{east}}, B{C{north}} or B{C{up}}.
        '''
        XyzLocal.__init__(self, enu, north, up, ltp=ltp, name=name)

#   @Property_RO
#   def xyzLocal(self):
#        return XyzLocal(self)


class Enu4Tuple(_NamedTuple):
    '''4-Tuple C{(east, north, up, ltp)}, in C{meter} except C{ltp}.
    '''
    _Names_ = (_east_, _north_, _up_,   _ltp_)
    _Units_ = ( Meter,  Meter,   Meter, _Pass)

    def _toEnu(self, Cls, Cls_kwds):
        '''(INTERNAL) Return C{Cls(..., **Cls_kwds)} instance.
        '''
        if issubclassof(Cls, XyzLocal):
            return Cls(*self, **_xkwds(Cls_kwds, name=self.name))
        else:
            return _4Tuple2Cls(self, Cls, Cls_kwds)

    @Property_RO
    def xyzLocal(self):
        '''Get this L{Enu4Tuple} as an L{XyzLocal}.
        '''
        return XyzLocal(*self, name=self.name)


class Local9Tuple(_NamedTuple):
    '''9-Tuple C{(x, y, z, lat, lon, height, ltp, ecef, M)} with I{local} C{x},
       C{y}, C{z} all in C{meter}, I{geodetic} C{lat}, C{lon}, C{height}, I{local
       tangent plane} C{ltp} (L{Ltp}), C{ecef} (L{Ecef9Tuple}) with I{geocentric}
       C{x}, C{y}, C{z}, I{geodetic} C{lat}, C{lon}, C{height} and I{concatenated}
       rotation matrix C{M} (L{EcefMatrix}) or C{None}.
    '''
    _Names_ = (_x_,   _y_,   _z_,   _lat_, _lon_, _height_, _ltp_, _ecef_, _M_)
    _Units_ = ( Meter, Meter, Meter, Lat,   Lon,   Height,  _Pass, _Pass,  _Pass)

    @Property_RO
    def azimuth(self):
        '''Get the I{local} Azimuth, bearing from North (C{degrees360}).
        '''
        return self.xyzLocal.aer4.azimuth

    @Property_RO
    def down(self):
        '''Get the I{local} Down, C{-z} component (C{meter}).
        '''
        return -self.z

    @Property_RO
    def east(self):
        '''Get the I{local} East, C{x} component (C{meter}).
        '''
        return self.x

    @Property_RO
    def elevation(self):
        '''Get the I{local} Elevation, tilt I{above} horizon (C{degrees90}).
        '''
        return self.xyzLocal.aer4.elevation

    @Property_RO
    def groundrange(self):
        '''Get the I{local} ground range, distance (C{meter}).
        '''
        return self.xyzLocal.aer4.groundrange

    @Property_RO
    def lam(self):
        '''Get the I{geodetic} longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @Property_RO
    def latlon(self):
        '''Get the I{geodetic} lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat, self.lon, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get the I{geodetic} lat-, longitude in C{degrees} and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

    @Property_RO
    def north(self):
        '''Get the I{local} North, C{y} component (C{meter}).
        '''
        return self.y

    @Property_RO
    def phi(self):
        '''Get the I{geodetic} latitude in C{radians} (C{float}).
        '''
        return self.philam.phi

    @Property_RO
    def philam(self):
        '''Get the I{geodetic} lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(radians(self.lat), radians(self.lon), name=self.name)

    @Property_RO
    def philamheight(self):
        '''Get the I{geodetic} lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    @Property_RO
    def slantrange(self):
        '''Get the I{local} slant Range, distance (C{meter}).
        '''
        return self.xyzLocal.aer4.slantrange

    def toAer(self, Aer=None, **Aer_kwds):
        '''Get the I{local} I{Azimuth, Elevation, slant Range} (AER) components.

           @kwarg Aer: Class to return AER (L{Aer}) or C{None}.
           @kwarg Aer_kwds: Optional, additional B{L{Aer}} keyword
                            arguments, ignored if B{C{Aer}} is C{None}.

           @return: AER as an L{Aer} instance or if C{B{Aer} is None},
                    an L{Aer4Tuple}C{(azimuth, elevation, slantrange, ltp)}.
        '''
        return self.xyzLocal.toAer(Aer=Aer, **Aer_kwds)

    def toCartesian(self, Cartesian=None, **Cartesian_kwds):
        '''Convert this I{local} to I{geocentric} C{(x, y, z)} (ECEF).

           @kwarg Cartesian: Optional class to return C{(x, y, z)} (C{Cartesian})
                             or C{None}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}} keyword
                                  arguments, ignored if C{B{Cartesian} is None}.

           @return: A C{B{Cartesian}(x, y, z, **B{Cartesian_kwds})} instance
                    or a L{Vector4Tuple}C{(x, y, z, h)} if C{B{Cartesian} is None}.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}
                             argument.
        '''
        return self.ecef.toCartesian(Cartesian=Cartesian, **Cartesian_kwds)  # PYCHOK _Tuple

    def toEnu(self, Enu=None, **Enu_kwds):
        '''Get the I{local} I{East, North, Up} (ENU) components.

           @kwarg Enu: Class to return ENU (L{Enu}) or C{None}.
           @kwarg Enu_kwds: Optional, additional B{L{Enu}} keyword
                            arguments, ignored if C{B{Enu} is None}.

           @return: ENU as an L{Enu} instance or if C{B{Enu} is None},
                    an L{Enu4Tuple}C{(east, north, up, ltp)}.
        '''
        return self.xyzLocal.toEnu(Enu=Enu, **Enu_kwds)

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Convert this I{local} to I{geodetic} C{(lat, lon, height)}.

           @kwarg LatLon: Optional class to return C{(lat, lon, height)}
                          (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: An instance of C{B{LatLon}(lat, lon, **B{LatLon_kwds})}
                    or if C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} respectively L{LatLon4Tuple}C{(lat, lon, height,
                    datum)} depending on whether C{datum} is un-/specified.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}
                             argument.
        '''
        return self.ecef.toLatLon(LatLon=LatLon, **LatLon_kwds)  # PYCHOK _Tuple

    def toNed(self, Ned=None, **Ned_kwds):
        '''Get the I{local} I{North, East, Down} (NED) components.

           @kwarg Ned: Class to return NED (L{Ned}) or C{None}.
           @kwarg Ned_kwds: Optional, additional B{L{Ned}} keyword
                            arguments, ignored if B{C{Ned}} is C{None}.

           @return: NED as an L{Ned} instance or if C{B{Ned} is None},
                    an L{Ned4Tuple}C{(north, east, down, ltp)}.
        '''
        return self.xyzLocal.toNed(Ned=Ned, **Ned_kwds)

    def toXyz(self, Xyz=None, **Xyz_kwds):
        '''Get the I{local} I{X, Y, Z} (XYZ) components.

           @kwarg Xyz: Class to return XYZ (L{XyzLocal}) or C{None}.
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: XYZ as an B{C{Xyz}} instance or if C{B{Xyz} is None},
                    an L{Xyz4Tuple}C{(x, y, z, ltp)}.
        '''
        return self.xyzLocal.toXyz(Xyz=Xyz, **Xyz_kwds)

    @Property_RO
    def up(self):
        '''Get the I{local} Up, C{z} component (C{meter}).
        '''
        return self.z

    @Property_RO
    def xyz(self):
        '''Get the I{local} C{(X, Y, Z)} components (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)  # like Ecef9Tuple.xyz

    @Property_RO
    def xyzLocal(self):
        '''Get this L{Local9Tuple} as an L{XyzLocal}.
        '''
        return XyzLocal(self)


_XyzLocals4 =  XyzLocal, Enu, Ned, Aer  # PYCHOK in .ltp
_XyzLocals5 = _XyzLocals4 + (Local9Tuple,)  # PYCHOK in .ltp


class Footprint5Tuple(_NamedTuple):
    '''5-Tuple C{(center, upperleft, upperight, loweright, lowerleft)}
       with the C{center} and 4 corners of a I{local} projection of a
       C{Frustum}, each an L{Xyz4Tuple}.

       @note: Misspelling of C{upperight} and C{loweright} is I{intentional}.
    '''
    _Names_ = ('center', 'upperleft', 'upperight', 'loweright', 'lowerleft')
    _Units_ = (_Pass,    _Pass,       _Pass,       _Pass,       _Pass)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
