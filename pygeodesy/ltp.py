
# -*- coding: utf-8 -*-

u'''I{Local cartesian} and I{local tangent plane} classes L{LocalCartesian} and L{Ltp},
L{LocalError} and L{Frustum}.

@see: U{Local tangent plane coordinates<https://WikiPedia.org/wiki/Local_tangent_plane_coordinates>}
      and class L{LocalCartesian}, transcribed from I{Charles Karney}'s C++ classU{LocalCartesian
      <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1LocalCartesian.html>}.
'''

# from pygeodesy.basics import issubclassof  # from ecef
from pygeodesy.datums import _WGS84, _xinstanceof
from pygeodesy.ecef import _EcefBase, EcefKarney, issubclassof, \
                           _llhn4, _xyzn4
from pygeodesy.errors import _TypesError, _ValueError
from pygeodesy.interns import EPS, NN, _ltp_, _M_, _lat0_, \
                             _lon0_, _name_, _0_, _0_0, _0_5, \
                             _2_0, _90_0, _180_0, _360_0
from pygeodesy.interns import _ecef_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.ltpTuples import Footprint5Tuple, Local9Tuple, \
                               _XyzLocals4, _XyzLocals5, Xyz4Tuple
from pygeodesy.named import _NamedBase
from pygeodesy.props import Property, Property_RO
from pygeodesy.units import Degrees, Meter
from pygeodesy.utily import sincos2d

from math import radians, tan

__all__ = _ALL_LAZY.ltp
__version__ = '21.04.17'

_Xyz_ = 'Xyz'


class LocalError(_ValueError):
    '''A L{LocalCartesian} or L{Ltp} related issue.
    '''
    pass


class LocalCartesian(_NamedBase):
    '''Conversion between geodetic C{(lat, lon, height)} and I{local cartesian}
       C{(x, y, z)} coordinates with I{geodetic} origin C{(lat0, lon0, height0)},
       transcribed from I{Karney}'s C++ class U{LocalCartesian
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1LocalCartesian.html>}.

       The C{z} axis is normal to the ellipsoid, the C{y} axis points due
       North.  The plane C{z = -heighth0} is tangent to the ellipsoid.

       The conversions all take place via geocentric coordinates using a
       geocentric L{EcefKarney}, by default the WGS84 datum/ellipsoid.

       @see: Class L{Ltp}.
    '''
    _ecef = EcefKarney(_WGS84)
    _t0   = None  # origin (..., lat0, lon0, height0, ...) L{Ecef9Tuple}

    def __init__(self, latlonh0=0, lon0=0, height0=0, ecef=None, name=NN):
        '''New L{LocalCartesian} converter.

           @kwarg latlonh0: Either a C{LatLon}, L{Ltp}, L{Ecef9Tuple} or
                            C{scalar} latitude of the (goedetic) origin
                            (C{degrees}).
           @kwarg lon0: Optional C{scalar} longitude of the (goedetic) origin
                        C{scalar} B{C{latlonh0}} (C{degrees}).
           @kwarg height0: Optional origin height (C{meter}), vertically
                           above (or below) the surface of the ellipsoid.
           @kwarg ecef: An ECEF converter (L{EcefKarney}).
           @kwarg name: Optional name (C{str}).

           @raise LocalError: If B{C{latlonh0}} not C{LatLon}, L{Ecef9Tuple},
                              C{scalar} or invalid or if B{C{lon0}} not
                              C{scalar} for C{scalar} B{C{latlonh0}} or
                              invalid or if B{C{height0}} invalid.

           @raise TypeError: Invalid B{C{ecef}}, not L{EcefKarney}.

           @note: If BC{latlonh0} is L{Ltp}, only the lat-, longitude and
                  height are duplicated, I{not} the ECEF converter.
        '''
        if ecef:
            _xinstanceof(EcefKarney, ecef=ecef)
            self._ecef = ecef
        self.reset(latlonh0, lon0, height0, name=name)

    @Property_RO
    def datum(self):
        '''Get the ECEF converter's datum (L{Datum}).
        '''
        return self.ecef.datum

    @Property_RO
    def ecef(self):
        '''Get the ECEF converter (L{EcefKarney}).
        '''
        return self._ecef

    def _ecef2local(self, ecef, Xyz, Xyz_kwds):
        '''(INTERNAL) Convert geocentric/geodetic to local, like I{forward}.

           @arg ecef: Geocentric (and geodetic) (L{Ecef9Tuple}).
           @arg Xyz: An L{XyzLocal}, L{Enu} or L{Ned} I{class} or C{None}.
           @arg Xyz_kwds: B{C{Xyz}} keyword arguments, ignored if C{B{Xyz}=None}.

           @return: An C{B{Xyz}(x, y, z, ltp, **B{Xyz_kwds}} instance or if
                    C{B{Xyz}=None}, an L{Local9Tuple}C{(x, y, z, lat, lon,
                    height, ltp, ecef, M)} with this C{ltp}, B{C{ecef}}
                    (L{Ecef9Tuple}) converted to this C{datum} and C{M=None},
                    always.
        '''
        ltp = self
        if ecef.datum != ltp.datum:
            ecef = ecef.toDatum(ltp.datum)
        x, y, z = self.M.rotate(ecef.xyz, *ltp._xyz0)
        r = Local9Tuple(x, y, z, ecef.lat, ecef.lon, ecef.height,
                                 ltp, ecef, None, name=ecef.name)
        if Xyz:
            if not issubclassof(Xyz, *_XyzLocals4):  # Vector3d
                raise _TypesError(_Xyz_, Xyz, *_XyzLocals4)
            r = r.toXyz(Xyz=Xyz, **Xyz_kwds)
        return r

    def forward(self, latlonh, lon=None, height=0, M=False):
        '''Convert I{geodetic} C{(lat, lon, height)} to I{local} cartesian
           C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, a L{Ltp}, an L{Ecef9Tuple} or
                         C{scalar} (geodetic) latitude (C{degrees}).
           @kwarg lon: Optional C{scalar} (geodetic) longitude for C{scalar}
                       B{C{latlonh}} (C{degrees}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: A L{Local9Tuple}C{(x, y, z, lat, lon, height, ltp, ecef, M)}
                    with I{local} C{x}, C{y}, C{z}, I{geodetic} C{(lat}, C{lon},
                    C{height}, this C{ltp}, C{ecef} (L{Ecef9Tuple}) with
                    I{geocentric} C{x}, C{y}, C{z} (and I{geodetic} C{lat},
                    C{lon}, C{height}) and the I{concatenated} rotation matrix
                    C{M} (L{EcefMatrix}) if requested.

           @raise LocalError: If B{C{latlonh}} not C{scalar}, C{LatLon}, L{Ltp},
                              L{Ecef9Tuple} or invalid or if B{C{lon}} not
                              C{scalar} for C{scalar} B{C{latlonh}} or invalid
                              or if B{C{height}} invalid.

           @see: Note at method L{EcefKarney.forward}.
        '''
        lat, lon, h, n = _llhn4(latlonh, lon, height, Error=LocalError)
        t = self.ecef.forward(lat, lon, h, M=M)
        x, y, z = self.M.rotate(t.xyz, *self._xyz0)
        m = self.M.multiply(t.M) if M else None
        return Local9Tuple(x, y, z, lat, lon, h, self, t, m, name=n or self.name)

    @Property_RO
    def height0(self):
        '''Get origin's height (C{meter}).
        '''
        return self._t0.height

    @Property_RO
    def lat0(self):
        '''Get origin's latitude (C{degrees}).
        '''
        return self._t0.lat

    def _local2ecef(self, local, nine=False, M=False):
        '''(INTERNAL) Convert I{local} to geocentric/geodetic, like I{.reverse}.

           @arg local: Local (L{XyzLocal}, L{Enu}, L{Ned}, L{Aer} or L{Local9Tuple}).
           @kwarg nine: Return 3- or 9-tuple (C{bool}).
           @kwarg M: Include the rotation matrix (C{bool}).

           @return: A I{geocentric} 3-tuple C{(x, y, z)} or if C{B{nine}=True},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)},
                    optionally including rotation matrix C{M} or C{None}.
        '''
        t = self.M.unrotate(local.xyz, *self._xyz0)
        if nine:
            t = self.ecef.reverse(*t, M=M)
        return t

    @Property_RO
    def lon0(self):
        '''Get origin's longitude (C{degrees}).
        '''
        return self._t0.lon

    @Property_RO
    def M(self):
        '''Get the rotation matrix (C{EcefMatrix}).
        '''
        return self._t0.M

    def reset(self, latlonh0=0, lon0=0, height0=0, name=NN):
        '''Reset the (geodetic) origin.

           @kwarg latlonh0: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                            latitude of the origin (C{degrees}).
           @kwarg lon0: Optional C{scalar} longitude of the origin for
                        C{scalar} B{C{latlonh0}} (C{degrees}).
           @kwarg height0: Optional origin height (C{meter}), vertically
                           above (or below) the surface of the ellipsoid.
           @kwarg name: Optional, new name (C{str}).

           @raise LocalError: If B{C{latlonh0}} not C{LatLon}, L{Ecef9Tuple},
                              C{scalar} or invalid or if B{C{lon0}} not
                              C{scalar} for C{scalar} B{C{latlonh0}} or
                              invalid or if B{C{height0}} invalid.
        '''
        self._update(True)  # force reset

        lat0, lon0, height0, n = _llhn4(latlonh0, lon0, height0,
                                        suffix=_0_, Error=LocalError)
        if name:
            n = name
        if n:
            self.rename(n)
            self.ecef.rename(n)
        self._t0 = self.ecef.forward(lat0, lon0, height0, M=True)

    def reverse(self, xyz, y=None, z=None, M=False):
        '''Convert I{local} C{(x, y, z)} to I{geodetic} C{(lat, lon, height)}.

           @arg xyz: A I{local} (L{XyzLocal}, L{Enu}, L{Ned}, L{Aer}, L{Local9Tuple}) or
                     local C{x} coordinate (C{scalar}).
           @kwarg y: Local C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: Local C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: Optionally, return the I{concatenated} rotation L{EcefMatrix} (C{bool}).

           @return: An L{Local9Tuple}C{(x, y, z, lat, lon, height, ltp, ecef, M)} with
                    I{local} C{x}, C{y}, C{z}, I{geodetic} C{lat}, C{lon}, C{height},
                    this C{ltp}, an C{ecef} (L{Ecef9Tuple}) with the I{geocentric} C{x},
                    C{y}, C{z} (and I{geodetic} C{lat}, C{lon}, C{height}) and the
                    I{concatenated} rotation matrix C{M} (L{EcefMatrix}) if requested.

           @raise LocalError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                              not C{scalar} for C{scalar} B{C{xyz}}.

           @see: Note at method L{EcefKarney.reverse}.
        '''
        x, y, z, n = _xyzn4(xyz, y, z, _XyzLocals5, Error=LocalError)
        xyz = self.M.unrotate((x, y, z), *self._xyz0)
        t = self.ecef.reverse(*xyz, M=M)
        m = self.M.multiply(t.M) if M else None
        return Local9Tuple(x, y, z, t.lat, t.lon, t.height, self, t, m, name=n or self.name)

    def toStr(self, prec=9):  # PYCHOK signature
        '''Return this L{LocalCartesian} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This L{LocalCartesian} representation (C{str}).
        '''
        return self.attrs(_lat0_, _lon0_, 'height0', _M_, 'ecef', _name_, prec=prec)

    @Property_RO
    def _xyz0(self):
        '''(INTERNAL) Get C{(x0, y0, z0)} as L{Vector3Tuple}.
        '''
        return self._t0.xyz


class Ltp(LocalCartesian):
    '''A I{local tangent plan} LTP, a sub-class of C{LocalCartesian} with
       configurable C{Ecef} and without optional rotation matrix.
    '''
    def __init__(self, latlonh0=0, lon0=0, height0=0, ecef=None, name=NN):
        '''New C{Ltp}.

           @kwarg latlonh0: Either a C{LatLon}, L{Ltp}, L{Ecef9Tuple} or
                            C{scalar} latitude of the (goedetic) origin
                            (C{degrees}).
           @kwarg lon0: Optional C{scalar} longitude of the (goedetic) origin
                        for C{scalar} B{C{latlonh0}} (C{degrees}).
           @kwarg height0: Optional origin height (C{meter}), vertically
                           above (or below) the surface of the ellipsoid.
           @kwarg ecef: Optional ECEF converter (L{EcefKarney}, l{EcefFarrell21},
                        L{EcefFarrell22}, L{EcefSudano}, L{EcefVeness} or
                        L{EcefYou} I{instance}), overriding default
                        L{EcefKarney}C{(datum=Datums.WGS84)}.
           @kwarg name: Optional name (C{str}).

           @return: New instance (C{Ltp}).

           @raise LocalError: If B{C{latlonh0}} not C{LatLon}, L{Ecef9Tuple},
                              C{scalar} or invalid or if B{C{lon0}} not
                              C{scalar} for C{scalar} B{C{latlonh0}} or
                              invalid or if B{C{height0}} invalid.

           @raise TypeError: Invalid B{C{ecef}}.
        '''
        if isinstance(latlonh0, Ltp):
            self._ecef = latlonh0.ecef
            self._t0   = latlonh0._to
            self.name  = name or latlonh0.name
        else:
            LocalCartesian.__init__(self, latlonh0, lon0=lon0, height0=height0, name=name)
        if ecef:
            self.ecef = ecef

    @Property
    def ecef(self):
        '''Get this LTP's ECEF converter (C{Ecef...} I{instance}).
        '''
        return self._ecef

    @ecef.setter  # PYCHOK setter!
    def ecef(self, ecef):
        '''Set this LTP's ECEF converter.

           @arg ecef: New ECEF converter (C{Ecef...} I{instance}).

           @raise TypeError: Invalid B{C{ecef}}.
        '''
        _xinstanceof(_EcefBase, ecef=ecef)
        if ecef != self._ecef:
            self._ecef = ecef
            self.reset(self._t0)


class Frustum(_NamedBase):
    '''A rectangular pyramid, typically representing a camera's I{field of view}
       (fov) and the intersection with (or projection to) a I{local tangent plane}.

       @see: U{Viewing frustum<https://WikiPedia.org/wiki/Viewing_frustum>}.
    '''
    _h_2     = _0_0   # half hfov
    _ltp     =  None  # local tangent plane
    _tan_h_2 = _0_0   # tan(_h_2)
    _v_2     = _0_0   # half vfov

    def __init__(self, hfov, vfov, ltp=None):
        '''New L{Frustum}.

           @arg hfov: Horizontal field of view (C{degrees180}).
           @arg vfov: Vertical field of view (C{degrees180}).
           @kwarg ltp: Optional I{local tangent plane} (L{Ltp}).

           @raise UnitError: Invalid B{C{hfov}} or B{C{vfov}}.

           @raise ValueError: Invalid B{C{hfov}} or B{C{vfov}}.
        '''
        self._h_2 = h = Degrees(hfov=hfov) * _0_5
        if not EPS < h < _90_0:
            raise _ValueError(hfov=hfov)

        self._v_2 = v = Degrees(vfov=vfov) * _0_5
        if not EPS < v < _90_0:
            raise _ValueError(vfov=vfov)

        self._tan_h_2 = t = tan(radians(h))
        if t < EPS:
            raise _ValueError(hfov=hfov)

        if ltp:
            self._ltp = _xLtp(ltp)

    @Property_RO
    def hfov(self):
        '''Get the horizontal C{fov} (C{degrees}).
        '''
        return Degrees(hfov=self._h_2 * _2_0)

    def footprint5(self, altitude, tilt, yaw=0, roll=0, z=_0_0, ltp=None):  # MCCABE 15
        '''Compute the center and corners of the intersection with (or projection
           to) the I{local tangent plane} (LTP).

           @arg altitude: Altitude (C{meter}) above I{local tangent plane}.
           @arg tilt: Pitch, elevation from horizontal (C{degrees180}), negative down.
           @kwarg yaw: Heading, bearing from North (C{degrees360}), clockwise from north.
           @kwarg roll: Roll, bank (C{degrees}), positive to the right and down.
           kwarg z: Optional height of the footprint (C{meter}) above I{local tangent plane}.
           @kwarg ltp: The I{local tangent plane} (L{Ltp}), overriding this
                       frustum's C{ltp}.

           @return: A L{Footprint5Tuple}C{(center, upperleft, upperight, loweright,
                    lowerleft)} with the C{center} and 4 corners each an L{Xyz4Tuple}.

           @raise TypeError: Invalid B{C{ltp}}.

           @raise UnitError: Invalid B{C{altitude}}, B{C{tilt}}, B{C{roll}} or B{C{z}}.

           @raise ValueError: If B{C{altitude}} too low or B{C{tilt}} or B{C{roll}}
                              -including B{C{vfow}} respectively B{C{hfov}}- over
                              the horizon.

           @see: U{Principal axes<https://WikiPedia.org/wiki/Aircraft_principal_axes>}.
        '''
        def _xy2(a, t, h_2, r, tan_h_2):
            # left and right corners, or swapped
            if r < EPS:  # no roll
                r = a * tan_h_2
                l = -r  # PYCHOK l is ell
            else:  # roll
                sl, cl, sr, cr = sincos2d(r + h_2, r - h_2)
                if cl < EPS or cr < EPS:
                    raise _ValueError(roll_hfov=r)
                r = -a * sr / cr  # negate right positive
                l = -a * sl / cl  # PYCHOK l is ell
            s, c = sincos2d(t)
            if abs(s) < EPS:
                raise _ValueError(tilt_vfov=t)
            y = a * c / s
            return (l, y), (r, y)

        def _xys(b, *xys):
            # rotate (x, y)'s clockwise
            s, c = sincos2d(b)
            for x, y in xys:
                yield (x * c + y * s), (y * c - x * s)

        a = Meter(altitude=altitude)
        if a < EPS:  # too low
            raise _ValueError(altitude=altitude)
        if z:
            z  = Meter(z=z)
            a -= z
            if a < EPS:  # z above a
                raise _ValueError(altitude_z=a)
        else:
            z = _0_0

        b =  Degrees(yaw=yaw) % _360_0
        t = -Degrees(tilt=tilt)
        if not EPS < t < _180_0:
            raise _ValueError(tilt=tilt)
        if t > _90_0:
            t =  _180_0 - t
            b = (_180_0 + b) % _360_0

        r = Degrees(roll=roll) % _360_0  # center
        if r:  # roll center
            s, c = sincos2d(r)
            if c < EPS:
                raise _ValueError(roll=r)
            x = -a * s / c
        else:
            x = _0_0
        # ground range
        s, c = sincos2d(t)
        if s < EPS:
            raise _ValueError(tilt=tilt)
        y = (a * c / s) if c > EPS else _0_0

        # center and corners, clockwise from upperleft
        xys = ((x, y),) + _xy2(a, t - self._v_2,  self._h_2, r,  self._tan_h_2) \
                        + _xy2(a, t + self._v_2, -self._h_2, r, -self._tan_h_2)  # swapped
        # rotate center and corners by yaw
        p = self.ltp if ltp is None else _xLtp(ltp)
        return Footprint5Tuple(*(Xyz4Tuple(x, y, z, p) for
                                           x, y in _xys(b, *xys)))

    @Property_RO
    def ltp(self):
        '''Get the I{local tangent plane} (L{Ltp}) or C{None}.
        '''
        return self._ltp

    @Property_RO
    def vfov(self):
        '''Get the vertical C{fov} (C{degrees}).
        '''
        return Degrees(vfov=self._v_2 * _2_0)


def _xLtp(ltp):
    '''(INTERNAL) Validate B{C{ltp}}.
    '''
    if isinstance(ltp, (LocalCartesian, Ltp)):
        return ltp
    raise _TypesError(_ltp_, ltp, Ltp, LocalCartesian)

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
