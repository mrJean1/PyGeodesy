
# -*- coding: utf-8 -*-

u'''I{Local Tangent Plane} (LTP) and I{local} cartesian coordinates.

I{Local cartesian} and I{local tangent plane} classes L{LocalCartesian} and L{Ltp},
L{LocalError} and L{Attitude} and L{Frustum}.

@see: U{Local tangent plane coordinates<https://WikiPedia.org/wiki/Local_tangent_plane_coordinates>}
      and class L{LocalCartesian}, transcoded from I{Charles Karney}'s C++ classU{LocalCartesian
      <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1LocalCartesian.html>}.
'''

from pygeodesy.basics import isscalar, issubclassof, map1, _xargs_names
from pygeodesy.constants import EPS, INT0, _umod_360, _0_0, _0_01, _0_0001, _0_5, \
                               _2_0, _60_0, _90_0, _100_0, _180_0, _3600_0, \
                               _N_1_0  # PYCHOK used!
from pygeodesy.datums import _WGS84, _xinstanceof
from pygeodesy.ecef import _EcefBase, EcefKarney, _llhn4, _xyzn4
from pygeodesy.errors import _NotImplementedError, _TypesError, _ValueError, _xkwds
from pygeodesy.fmath import fdot, Fhorner
from pygeodesy.fsums import _floor, Fsum, fsum_, fsum1_
from pygeodesy.interns import NN, _0_, _COMMASPACE_, _DOT_, _ecef_, _height_, \
                             _invalid_, _lat0_, _lon0_, _ltp_, _M_, _name_, _too_
# from pygeodesy.lazily import _ALL_LAZY  # from .vector3d
from pygeodesy.ltpTuples import Attitude4Tuple, ChLV9Tuple, Footprint5Tuple, \
                                Local9Tuple, _XyzLocals4, _XyzLocals5, Xyz4Tuple
from pygeodesy.named import _NamedBase, notOverloaded
from pygeodesy.namedTuples import LatLon3Tuple, LatLon4Tuple, Vector3Tuple
from pygeodesy.props import Property, Property_RO, property_doc_, property_RO, \
                           _update_all
from pygeodesy.streprs import Fmt, strs, unstr
from pygeodesy.units import Bearing, Degrees, Meter
from pygeodesy.utily import cotd, sincos2d, sincos2d_, tand, tand_, wrap180, wrap360
from pygeodesy.vector3d import _ALL_LAZY, Vector3d

# from math import floor as _floor  # from .fsums

__all__ = _ALL_LAZY.ltp
__version__ = '22.09.29'

_height0_ = _height_ + _0_
_narrow_  = 'narrow'
_wide_    = 'wide'
_Xyz_     = 'Xyz'


def _fov_2(**fov):
    # Half a field-of-view angle in C{degrees}.
    f = Degrees(Error=LocalError, **fov) * _0_5
    if EPS < f < _90_0:
        return f
    t = _invalid_ if f < 0 else _too_(_wide_ if f > EPS else _narrow_)
    raise LocalError(txt=t, **fov)


class Attitude(_NamedBase):
    '''The orientation of a plane or camera in space.
    '''
    _alt  = Meter(  alt =_0_0)
    _roll = Degrees(roll=_0_0)
    _tilt = Degrees(tilt=_0_0)
    _yaw  = Bearing(yaw =_0_0)

    def __init__(self, alt_attitude=INT0, tilt=INT0, yaw=INT0, roll=INT0, name=NN):
        '''New L{Attitude}.

           @kwarg alt_attitude: An altitude (C{meter}) above earth or an attitude
                                (L{Attitude} or L{Attitude4Tuple}) with the
                                C{B{alt}itude}, B{C{tilt}}, B{C{yaw}} and B{C{roll}}.
           @kwarg tilt: Pitch, elevation from horizontal (C{degrees180}), negative down
                        (clockwise rotation along and around the x- or East axis).
           @kwarg yaw: Bearing, heading (compass C{degrees360}), clockwise from North
                       (counter-clockwise rotation along and around the z- or Up axis).
           @kwarg roll: Roll, bank (C{degrees180}), positive to the right and down
                        (clockwise rotation along and around the y- or North axis).
           @kwarg name: Optional name C{str}).

           @raise AttitudeError: Invalid B{C{alt_attitude}}, B{C{tilt}}, B{C{yaw}} or
                                 B{C{roll}}.

           @see: U{Principal axes<https://WikiPedia.org/wiki/Aircraft_principal_axes>} and
                 U{Yaw, pitch, and roll rotations<http://Planning.CS.UIUC.edu/node102.html>}.
        '''
        if isscalar(alt_attitude):
            t = Attitude4Tuple(alt_attitude, tilt, yaw, roll)
        else:
            try:
                t = alt_attitude.atyr
            except AttributeError:
                raise AttitudeError(alt=alt_attitude, tilt=tilt, yaw=yaw, rol=roll)
        for n, v in t.items():
            if v:
                setattr(self, n, v)
        n = name or t.name
        if n:
            self.name = n

    @property_doc_(' altitude above earth in C{meter}.')
    def alt(self):
        return self._alt

    @alt.setter  # PYCHOK setter!
    def alt(self, alt):  # PYCHOK no cover
        a = Meter(alt=alt, Error=AttitudeError)
        if self._alt != a:
            _update_all(self)
            self._alt = a

    altitude = alt

    @Property_RO
    def atyr(self):
        '''Return this attitude's alt[itude], tilt, yaw and roll as an L{Attitude4Tuple}.
        '''
        return Attitude4Tuple(self.alt, self.tilt, self.yaw, self.roll, name=self.name)

    @Property_RO
    def matrix(self):
        '''Get the 3x3 rotation matrix C{R(yaw)·R(tilt)·R(roll)}, aka I{ZYX} (C{float}, row-order).

           @see: The matrix M of case 10 in U{Appendix A
                 <https://ntrs.NASA.gov/api/citations/19770019231/downloads/19770019231.pdf>}.
        '''
        def _5to3(x, y, _y, z, _z):
            return x, fsum1_(y, _y), fsum1_(z, _z)

        r0, r1, r2 = self._rows3
        return _5to3(*r0), _5to3(*r1), r2

    @property_doc_(' roll/bank in C{degrees180}, positive to the right and down.')
    def roll(self):
        return self._roll

    @roll.setter  # PYCHOK setter!
    def roll(self, roll):
        r = Degrees(roll=roll, wrap=wrap180, Error=AttitudeError)
        if self._roll != r:
            _update_all(self)
            self._roll = r

    bank = roll

    @Property_RO
    def _rows3(self):
        # to follow the definitions of rotation angles alpha, beta and gamma:
        # negate yaw since yaw is counter-clockwise around the z-axis, swap
        # tilt and roll since tilt is around the x- and roll around the y-axis
        sa, ca, sb, cb, sg, cg = sincos2d_(-self.yaw, self.roll, self.tilt)
        return ((ca * cb,  ca * sb * sg, -sa * cg,  ca * sb * cg,  sa * sg),
                (sa * cb,  sa * sb * sg,  ca * cg,  sa * sb * cg, -ca * sg),
                (    -sb,       cb * sg,                 cb * cg))

    def rotate(self, x_xyz, y=None, z=None, Vector=None, **Vector_kwds):
        '''Transform a (local) cartesian by this attitude's matrix.

           @arg x_xyz: X component of vector (C{scalar}) or (3-D) vector
                       (C{Cartesian}, L{Vector3d} or L{Vector3Tuple}).
           @kwarg y: Y component of vector (C{scalar}), same units as B{C{x}}.
           @kwarg z: Z component of vector (C{scalar}), same units as B{C{x}}.
           @kwarg Vector: Class to return transformed point (C{Cartesian},
                          L{Vector3d} or C{Vector3Tuple}) or C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                               ignored if C{B{Vector} is None}.

           @return: A B{C{Vector}} instance or a L{Vector3Tuple}C{(x, y, z)} if
                    C{B{Vector}=None}.

           @see: U{Yaw, pitch, and roll rotations<http://Planning.CS.UIUC.edu/node102.html>}.
        '''
        try:
            x, y, z = map( float, x_xyz.xyz)
        except AttributeError:
            x, y, z = map1(float, x_xyz, y, z)

        r0, r1, r2 = self._rows3
        X = fdot(r0, x, y, y, z, z)
        Y = fdot(r1, x, y, y, z, z)
        Z = fdot(r2, x, y, z)
        return Vector3Tuple(X, Y, Z, name=self.name) if Vector is None else \
                     Vector(X, Y, Z, **_xkwds(Vector_kwds, name=self.name))

    @property_doc_(' tilt/pitch/elevation from horizontal in C{degrees180}, negative down.')
    def tilt(self):
        return self._tilt

    @tilt.setter  # PYCHOK setter!
    def tilt(self, tilt):
        t = Degrees(tilt=tilt, wrap=wrap180, Error=AttitudeError)
        if self._tilt != t:
            _update_all(self)
            self._tilt = t

    elevation = pitch = tilt

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Format this attitude as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: This attitude (C{str}).
        '''
        return self.atyr.toStr(prec=prec, sep=sep)

    @Property_RO
    def tyr3d(self):
        '''Get this attitude's (3-D) directional vector (L{Vector3d}).

           @see: U{Yaw, pitch, and roll rotations<http://Planning.CS.UIUC.edu/node102.html>}.
        '''
        def _r2d(r):
            return fsum_(_N_1_0, *r)

        return Vector3d(*map1(_r2d, *self._rows3), name=tyr3d.__name__)

    @property_doc_(' yaw/bearing/heading in compass C{degrees360}, clockwise from North.')
    def yaw(self):
        return self._yaw

    @yaw.setter  # PYCHOK setter!
    def yaw(self, yaw):
        y = Bearing(yaw=yaw, Error=AttitudeError)
        if self._yaw != y:
            _update_all(self)
            self._yaw = y

    bearing = heading = yaw


class AttitudeError(_ValueError):
    '''An L{Attitude} or L{Attitude4Tuple} issue.
    '''
    pass


class Frustum(_NamedBase):
    '''A rectangular pyramid, typically representing a camera's I{field-of-view}
       (fov) and the intersection with (or projection to) a I{local tangent plane}.

       @see: U{Viewing frustum<https://WikiPedia.org/wiki/Viewing_frustum>}.
    '''
    _h_2     = _0_0   # half hfov in degrees
    _ltp     =  None  # local tangent plane
    _tan_h_2 = _0_0   # tan(_h_2)
    _v_2     = _0_0   # half vfov in degrees

    def __init__(self, hfov, vfov, ltp=None):
        '''New L{Frustum}.

           @arg hfov: Horizontal field-of-view (C{degrees180}).
           @arg vfov: Vertical field-of-view (C{degrees180}).
           @kwarg ltp: Optional I{local tangent plane} (L{Ltp}).

           @raise LocalError: Invalid B{C{hfov}} or B{C{vfov}}.
        '''
        self._h_2 = h = _fov_2(hfov=hfov)
        self._v_2 =     _fov_2(vfov=vfov)

        self._tan_h_2 = tand(h, fov_2=h)

        if ltp:
            self._ltp = _xLtp(ltp)

    def footprint5(self, alt_attitude, tilt=0, yaw=0, roll=0, z=_0_0, ltp=None):  # MCCABE 15
        '''Compute the center and corners of the intersection with (or projection
           to) the I{local tangent plane} (LTP).

           @arg alt_attitude: An altitude (C{meter}) above I{local tangent plane} or
                              an attitude (L{Attitude} or L{Attitude4Tuple}) with the
                              C{B{alt}itude}, B{C{tilt}}, B{C{yaw}} and B{C{roll}}.
           @kwarg tilt: Pitch, elevation from horizontal (C{degrees}), negative down
                        (clockwise rotation along and around the x- or East axis).
           @kwarg yaw: Bearing, heading (compass C{degrees}), clockwise from North
                       (counter-clockwise rotation along and around the z- or Up axis).
           @kwarg roll: Roll, bank (C{degrees}), positive to the right and down
                        (clockwise rotation along and around the y- or North axis).
           @kwarg z: Optional height of the footprint (C{meter}) above I{local tangent plane}.
           @kwarg ltp: The I{local tangent plane} (L{Ltp}), overriding this
                       frustum's C{ltp}.

           @return: A L{Footprint5Tuple}C{(center, upperleft, upperight, loweright,
                    lowerleft)} with the C{center} and 4 corners, each an L{Xyz4Tuple}.

           @raise TypeError: Invalid B{C{ltp}}.

           @raise UnitError: Invalid B{C{altitude}}, B{C{tilt}}, B{C{roll}} or B{C{z}}.

           @raise ValueError: If B{C{altitude}} too low, B{C{z}} too high or B{C{tilt}}
                              or B{C{roll}} -including B{C{vfov}} respectively B{C{hfov}}-
                              over the horizon.

           @see: U{Principal axes<https://WikiPedia.org/wiki/Aircraft_principal_axes>}.
        '''
        def _xy2(a, e, h_2, tan_h_2, r):
            # left and right corners, or swapped
            if r < EPS:  # no roll
                r =  a * tan_h_2
                l = -r  # PYCHOK l is ell
            else:  # roll
                r, l = tand_(r - h_2, r + h_2, roll_hfov=r)  # PYCHOK l is ell
                r *= -a  # negate right positive
                l *= -a  # PYCHOK l is ell
            y = a * cotd(e, tilt_vfov=e)
            return (l, y), (r, y)

        def _xyz5(b, xy5, z, ltp):
            # rotate (x, y)'s by bearing, clockwise
            s, c = sincos2d(b)
            for x, y in xy5:
                yield Xyz4Tuple(fsum1_(x * c,  y * s),
                                fsum1_(y * c, -x * s), z, ltp)

        try:
            a, t, y, r = alt_attitude.atyr
        except AttributeError:
            a, t, y, r = alt_attitude, tilt, yaw, roll

        a = Meter(altitude=a)
        if a < EPS:  # too low
            raise _ValueError(altitude=a)
        if z:  # PYCHOK no cover
            z  = Meter(z=z)
            a -= z
            if a < EPS:  # z above a
                raise _ValueError(altitude_z=a)
        else:
            z = _0_0

        b =  Degrees(yaw=y, wrap=wrap360)  # bearing
        e = -Degrees(tilt=t, wrap=wrap180)  # elevation, pitch
        if not EPS < e < _180_0:
            raise _ValueError(tilt=t)
        if e > _90_0:
            e = _180_0 - e
            b = _umod_360(b + _180_0)

        r = Degrees(roll=r, wrap=wrap180)  # roll center
        x = (-a * tand(r, roll=r)) if r else _0_0
        y =   a * cotd(e, tilt=t)  # ground range
        if abs(y) < EPS:
            y = _0_0

        # center and corners, clockwise from upperleft, rolled
        xy5 = ((x, y),) + _xy2(a, e - self._v_2,  self._h_2,  self._tan_h_2, r) \
                        + _xy2(a, e + self._v_2, -self._h_2, -self._tan_h_2, r)  # swapped
        # turn center and corners by yaw, clockwise
        p = self.ltp if ltp is None else ltp  # None OK
        return Footprint5Tuple(_xyz5(b, xy5, z, p))  # *_xyz5

    @Property_RO
    def hfov(self):
        '''Get the horizontal C{fov} (C{degrees}).
        '''
        return Degrees(hfov=self._h_2 * _2_0)

    @Property_RO
    def ltp(self):
        '''Get the I{local tangent plane} (L{Ltp}) or C{None}.
        '''
        return self._ltp

    def toStr(self, prec=3, fmt=Fmt.F, sep=_COMMASPACE_):  # PYCHOK signature
        '''Convert this frustum to a "hfov, vfov, ltp" string.

           @kwarg prec: Number of (decimal) digits, unstripped (0..8 or C{None}).
           @kwarg fmt: Optional, C{float} format (C{str}).
           @kwarg sep: Separator to join (C{str}).

           @return: Frustum in the specified form (C{str}).
        '''
        t = self.hfov, self.vfov
        if self.ltp:
            t += self.ltp,
        t = strs(t, prec=prec, fmt=fmt)
        return sep.join(t) if sep else t

    @Property_RO
    def vfov(self):
        '''Get the vertical C{fov} (C{degrees}).
        '''
        return Degrees(vfov=self._v_2 * _2_0)


class LocalError(_ValueError):
    '''A L{LocalCartesian} or L{Ltp} related issue.
    '''
    pass


class LocalCartesian(_NamedBase):
    '''Conversion between geodetic C{(lat, lon, height)} and I{local cartesian}
       C{(x, y, z)} coordinates with I{geodetic} origin C{(lat0, lon0, height0)},
       transcoded from I{Karney}'s C++ class U{LocalCartesian
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1LocalCartesian.html>}.

       The C{z} axis is normal to the ellipsoid, the C{y} axis points due
       North.  The plane C{z = -height0} is tangent to the ellipsoid.

       The conversions all take place via geocentric coordinates using a
       geocentric L{EcefKarney}, by default the WGS84 datum/ellipsoid.

       @see: Class L{Ltp}.
    '''
    _ecef   = EcefKarney(_WGS84)
    _t0     = None  # origin (..., lat0, lon0, height0, ...) L{Ecef9Tuple}
    _9Tuple = Local9Tuple

    def __init__(self, latlonh0=INT0, lon0=INT0, height0=INT0, ecef=None, name=NN):
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
        if isinstance(latlonh0, LocalCartesian):
            self._ecef = latlonh0.ecef
            self._t0   = latlonh0._t0
            self.name  = name or latlonh0.name
        else:
            self.reset(latlonh0, lon0, height0, name=name)
        if ecef:  # PYCHOK no cover
            _xinstanceof(EcefKarney, ecef=ecef)
            self._ecef = ecef

    def __eq__(self, other):
        '''Compare this and an other instance.

           @arg other: The other ellipsoid (L{LocalCartesian} or L{Ltp}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return other is self or (isinstance(other, self.__class__) and
                                     other.ecef == self.ecef and
                                     other._t0  == self._t0)

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
           @arg Xyz_kwds: B{C{Xyz}} keyword arguments, ignored if C{B{Xyz} is None}.

           @return: An C{B{Xyz}(x, y, z, ltp, **B{Xyz_kwds}} instance or if
                    C{B{Xyz} is None}, an L{Local9Tuple}C{(x, y, z, lat, lon,
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

    def forward(self, latlonh, lon=None, height=0, M=False, name=NN):
        '''Convert I{geodetic} C{(lat, lon, height)} to I{local} cartesian
           C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, a L{Ltp}, an L{Ecef9Tuple} or
                         C{scalar} (geodetic) latitude (C{degrees}).
           @kwarg lon: Optional C{scalar} (geodetic) longitude for C{scalar}
                       B{C{latlonh}} (C{degrees}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional name (C{str}).

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
        '''
        lat, lon, h, n = _llhn4(latlonh, lon, height, Error=LocalError, name=name)
        t = self.ecef._forward(lat, lon, h, n, M=M)
        x, y, z = self.M.rotate(t.xyz, *self._xyz0)
        m = self.M.multiply(t.M) if M else None
        return self._9Tuple(x, y, z, lat, lon, h, self, t, m, name=n or self.name)

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

    @Property_RO
    def latlonheight0(self):
        '''Get the origin's lat-, longitude and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return LatLon3Tuple(self.lat0, self.lon0, self.height0, name=self.name)

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

    def reset(self, latlonh0=INT0, lon0=INT0, height0=INT0, name=NN):
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
        _update_all(self)  # force reset

        lat0, lon0, height0, n = _llhn4(latlonh0, lon0, height0,
                                        suffix=_0_, Error=LocalError, name=name)
        if n:
            self.rename(n)
        else:
            n = self.name
        self._t0 = self.ecef._forward(lat0, lon0, height0, n, M=True)

    def reverse(self, xyz, y=None, z=None, M=False, name=NN):
        '''Convert I{local} C{(x, y, z)} to I{geodetic} C{(lat, lon, height)}.

           @arg xyz: A I{local} (L{XyzLocal}, L{Enu}, L{Ned}, L{Aer}, L{Local9Tuple}) or
                     local C{x} coordinate (C{scalar}).
           @kwarg y: Local C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: Local C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: Optionally, return the I{concatenated} rotation L{EcefMatrix},
                     I{iff avaialble} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: An L{Local9Tuple}C{(x, y, z, lat, lon, height, ltp, ecef, M)} with
                    I{local} C{x}, C{y}, C{z}, I{geodetic} C{lat}, C{lon}, C{height},
                    this C{ltp}, an C{ecef} (L{Ecef9Tuple}) with the I{geocentric} C{x},
                    C{y}, C{z} (and I{geodetic} C{lat}, C{lon}, C{height}) and the
                    I{concatenated} rotation matrix C{M} (L{EcefMatrix}) if requested.

           @raise LocalError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                              not C{scalar} for C{scalar} B{C{xyz}}.
        '''
        x, y, z, n = _xyzn4(xyz, y, z, _XyzLocals5, Error=LocalError, name=name)
        c = self.M.unrotate((x, y, z), *self._xyz0)
        t = self.ecef.reverse(*c, M=M)
        m = self.M.multiply(t.M) if M else None
        return self._9Tuple(x, y, z, t.lat, t.lon, t.height, self, t, m, name=n or self.name)

    def toStr(self, prec=9, **unused):  # PYCHOK signature
        '''Return this L{LocalCartesian} as a string.

           @kwarg prec: Precision, number of (decimal) digits (0..9).

           @return: This L{LocalCartesian} representation (C{str}).
        '''
        return self.attrs(_lat0_, _lon0_, _height0_, _M_, _ecef_, _name_, prec=prec)

    @Property_RO
    def _xyz0(self):
        '''(INTERNAL) Get C{(x0, y0, z0)} as L{Vector3Tuple}.
        '''
        return self._t0.xyz


class Ltp(LocalCartesian):
    '''A I{local tangent plan} LTP, a sub-class of C{LocalCartesian} with
       configurable ECEF converter and without optional rotation matrix.
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
        '''Set this LTP's ECEF converter (C{Ecef...} I{instance}).

           @raise TypeError: Invalid B{C{ecef}}.
        '''
        _xinstanceof(_EcefBase, ecef=ecef)
        if ecef != self._ecef:  # PYCHOK no cover
            self.reset(self._t0)
            self._ecef = ecef


class _ChLV(object):
    '''(INTERNAL) Base class for C{ChLV*} classes.
    '''
    _YXm = 1.e-6  # Swiss a, b per Y, X meter

    def forward(self, latlonh, lon=None, height=0, M=None, name=NN):
        '''Convert WGS84 geodetic to I{Swiss} projection coordinates.

           @arg latlonh: Either a C{LatLon}, L{Ltp} or C{scalar} (geodetic) latitude
                         (C{degrees}).
           @kwarg lon: Optional, C{scalar} (geodetic) longitude for C{scalar} B{C{latlonh}}
                       (C{degrees}).
           @kwarg height: Optional, height, vertically above (or below) the surface of the
                          ellipsoid (C{meter}) for C{scalar} B{C{latlonh}}.
           @kwarg M: Optionally, return the I{concatenated} rotation L{EcefMatrix}, I{iff
                     avaialble} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: A L{ChLV9Tuple}C{(Y, X, h_, lat, lon, height, ltp, ecef, M)} with
                    I{unfalsed Swiss Y, X, h_} coordinates and height, the given I{geodetic}
                    C{(lat}, C{lon}, C{height}, this C{ChLV*} and C{ecef} (L{Ecef9Tuple})
                    I{at Bern, Ch} and rotation matrix C{M}.

           @note: The returned C{ltp} is this C{ChLV}, C{ChLVa} or C{ChLVe} instance.

           @raise LocalError: If B{C{latlonh}} not C{scalar}, C{LatLon} or L{Ltp},
                              invalid or if B{C{lon}} not C{scalar} for C{scalar}
                              B{C{latlonh}} or invalid or if B{C{height}} invalid.
        '''
        notOverloaded(self, latlonh, lon=lon, height=height, M=M, name=name)

    def reverse(self, enh_, n=None, h_=0, M=None, name=NN):
        '''Convert I{Swiss} projection to WGS84 geodetic coordinates.

           @arg enh_: A Swiss projection (L{ChLV9Tuple}) or the C{scalar}, I{falsed Swiss
                      E_LV95} or I{y_LV03} coordinate (C{meter}).
           @kwarg n: I{Falsed Swiss N_LV85} or I{x_LV03} coordinate for C{scalar}
                     B{C{enh_}} and B{C{h_}} (C{meter}).
           @kwarg h_: I{Swiss h'} height for C{scalar} B{C{enh_}} and B{C{n}} (C{meter}).
           @kwarg M: Optionally, return the I{concatenated} rotation L{EcefMatrix}, I{iff
                     avaialble} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: A L{ChLV9Tuple}C{(Y, X, h_, lat, lon, height, ltp, ecef, M)} with
                    the I{unfalsed Swiss (Y, X, h_)} coordinates and height, I{geodetic}
                    C{(lat, lon, height)}, this C{ChLV*} and C{ecef} (L{Ecef9Tuple})
                    I{at Bern, Ch} and with rotation matrix C{B{M}=None}.

           @note: The returned C{ltp} is this C{ChLV}, C{ChLVa} or C{ChLVe} instance.

           @raise LocalError: Invalid B{C{enh_}} or B{C{n}} or B{C{h_}} not C{scalar}
                              for C{scalar} B{C{enh_}.
        '''
        notOverloaded(self, enh_, n=n, h_=h_, M=M, name=name)

    @property_RO
    def _enh_n_h(self):
        '''(INTERNAL) Get C{ChLV*.reverse}'s args[1:4] names, I{once}.
        '''
        ChLV._enh_n_h = t = _xargs_names(ChLV.reverse)[1:4]  # overwrite property
        # assert _xargs_names(ChLVa.reverse)[1:4] == t
        # assert _xargs_names(ChLVe.reverse)[1:4] == t
        return t

    @staticmethod
    def _falsing2(Y, X, LV95, add):
        '''(INTERNAL) Return the e- and n-falsing.
        '''
        if LV95 is None:
            e = n = 0
        else:
            e = 2600000.0 if LV95 else 600000.0
            n = 1200000.0 if LV95 else 200000.0
        return (Y + e, X + n) if add else (Y - e, X - n)

    @staticmethod
    def _llh2abh_3(lat, lon, h):
        '''(INTERNAL) Helper for C{ChLVa/e.forward}.
        '''
        def _dms(ds, p, q, swap):
            d = _floor(ds)
            t = (ds - d) * p
            m = _floor(t)
            s = (t  - m) * p
            if swap:
                d, s = s, d
            return d + (m + s * q) * q

        def _deg2ab(deg, arcs):
            s = _dms(deg,  _60_0, _0_01, False)  # dec2sexag
            s = _dms(  s, _100_0, _60_0, True)   # sexag2sec
            return (s - arcs) * _0_0001

        a  = _deg2ab(lat, ChLV._sLat)  # phi', lat_aux
        b  = _deg2ab(lon, ChLV._sLon)  # lam', lng_aux
        h_ =  fsum_(h, -ChLV.Bern.height, 2.73 * b, 6.94 * a)
        return a, b, h_

    @staticmethod
    def _YXh_2abh3(Y, X, h_):
        '''(INTERNAL) Helper for C{ChLVa/e.reverse}.
        '''
        def _YX2ab(YX):
            return YX * _ChLV._YXm

        a, b = map1(_YX2ab, Y, X)
        h = fsum_(h_, ChLV.Bern.height, -12.6 * a, -22.64 * b)
        return a, b, h

    def _YXh_n4(self, enh_, n, h_, name):
        '''(INTERNAL) Helper for C{ChLV*.reverse}.
        '''
        Y, X, h_, name = _xyzn4(enh_, n, h_, ChLV9Tuple, name=name,
                                            _xyz_y_z_names=self._enh_n_h)
        if isinstance(enh_, ChLV9Tuple):
            Y, X = enh_.Y, enh_.X
        else:  # isscalar(enh_)
            Y, X = ChLV.unfalse2(Y, X)
        return Y, X, h_, name


class ChLV(_ChLV, Ltp):
    '''Conversion between I{WGS84 geodetic} and I{Swiss} projection coordinates
       using L{pygeodesy.EcefKarney}'s Earth-Centered, Earth-Fixed (ECEF) methods.

       @see: U{Swiss projection formulas<<https://www.SwissTopo.admin.CH/en/
             maps-data-online/calculation-services.html>}, pp 7-9.
    '''
    _9Tuple = ChLV9Tuple

    _sLat = 169028.66  # Bern, Ch in ...
    _sLon =  26782.5   # ... arc-seconds
    # lat, lon, height == 46°57'08.66", 7°26'22.50", 49.55m ("new" 46°57'07.89", 7°26'22.335")
    Bern  = LatLon4Tuple(_sLat / _3600_0, _sLon / _3600_0, 49.55, _WGS84, name='Bern')

    def __init__(self, latlonh0=Bern, **Ltp_kwds):
        '''New C{ChLV} L{pygeodesy.EcefKarney} I{WGS84-Swiss} converter, centered at I{Bern, Ch}.

           @see: L{Ltp.__init__} for more information.
        '''
        Ltp.__init__(self, latlonh0, **_xkwds(Ltp_kwds, ecef=None, name=ChLV.Bern.name))

    def forward(self, latlonh, lon=None, height=0, M=None, name=NN):  # PYCHOK unused M
        # overloaded for the _ChLV.forward.__doc__
        return Ltp.forward(self, latlonh, lon=lon, height=height, M=M, name=name)

    def reverse(self, enh_, n=None, h_=0, M=None, name=NN):  # PYCHOK unused M
        # overloaded for the _ChLV.reverse.__doc__
        Y, X, h_, name = self._YXh_n4(enh_, n, h_, name=name)
        return Ltp.reverse(self, Y, X, h_, M=M, name=name)

    @staticmethod
    def false2(Y, X, LV95=True):
        '''Add the C{LV95} or C{LV03} falsing.

           @arg Y: I{Unfalsed Swiss Y} coordinate, easting (C{meter}).
           @arg X: I{Unfalsed Swiss X} coordinate, northing (C{meter}).
           @kwarg LV95: Add C{LV95} or C{LV03} if C{True} respectively
                        C{False}, otherwise don't false.

           @return: 2-Tuple C{(e, n)} with falsed B{C{Y}} respectively B{C{X}}.
        '''
        return _ChLV._falsing2(Y, X, LV95, True)

    @staticmethod
    def isLV95(e, n, raiser=True):
        '''Determine the falsing as C{LV95}, C{LV03} or neither.

           @arg e: I{Falsed Swiss E_LV95} or I{y_LV03} coordinate (C{meter}).
           @arg n: I{Falsed Swiss N_LV95} or I{x_LV03} coordinate (C{meter}).
           @kwarg raiser: If C{True}, throw a L{LocalError} if B{C{e}} and
                          B{C{n}} are invalid C{LV95} and LC{LV03}.

           @return: C{True} or C{False} if B{C{e}} and B{C{n}} are valid C{LV95}
                    respectively C{LV03} coordinates, C{None} otherwise.
        '''
        # @see: U{Map<https://www.SwissTopo.admin.CH/en/knowledge-facts/
        #       surveying-geodesy/reference-frames/local/lv95.html>}
        def _isLV03(e, n):
            return 400000 < e < 900000 and 90000 < n < 400000

        if _isLV03(e, n):
            return False
        elif _isLV03(e - 2000000, n - 1000000):
            return True
        if raiser:
            t = unstr(ChLV.isLV95.__name__, e=e, n=n)
            raise LocalError(t)
        return None

    @staticmethod
    def unfalse2(e, n, LV95=None):
        '''Remove the C{LV95} or C{LV03} falsing.

           @arg e: I{Falsed Swiss E_LV95} or I{y_LV03} coordinate (C{meter}).
           @arg n: I{Falsed Swiss N_LV95} or I{x_LV03} coordinate (C{meter}).
           @kwarg LV95: Remove the C{LV95} falsing if C{True}, remove the
                        C{LV03} falsing if C{False}, otherwise use method
                        C{isLV95(B{e}, B{n})}.

           @return: 2-Tuple C{(Y, X)}, unfalsed B{C{e}}, respectively B{C{n}}.
        '''
        if LV95 is None:
            LV95 = ChLV.isLV95(e, n)
        return _ChLV._falsing2(e, n, LV95, False)


class ChLVa(_ChLV, LocalCartesian):
    '''Conversion between I{WGS84 geodetic} and I{Swiss} projection coordinates
       using the U{Approximate<https://www.SwissTopo.admin.CH/en/maps-data-online/
       calculation-services.html>} formulas, pp 13-14.

       @see: Older U{references<https://GitHub.com/alphasldiallo/Swisstopo-WGS84-LV03>}.
    '''
    def __init__(self, name=ChLV.Bern.name):
        '''New C{ChLV} I{Approximate WGS84-Swiss} converter, centered at I{Bern, Ch}.

           @kwarg name: Optional name (C{str}), overriding C{Bern.name}.
        '''
        LocalCartesian.__init__(self, latlonh0=ChLV.Bern, name=name)

    def forward(self, latlonh, lon=None, height=0, M=None, name=NN):
        # overloaded for the _ChLV.forward.__doc__
        lat, lon, h, name = _llhn4(latlonh, lon, height, name=name)

        a,  b, h_ = _ChLV._llh2abh_3(lat, lon, h)
        a2, b2    =  a**2, b**2

        Y = fsum_( 72.37, 211455.93 * b,
                          -10938.51 * b * a,
                              -0.36 * b * a2,
                             -44.54 * b * b2)  # + 600_000
        X = fsum_(147.07, 308807.95 * a,
                            3745.25 * b2,
                              76.63 * a2,
                            -194.56 * b2 * a,
                             119.79 * a2 * a)  # + 200_000
        return self._ChLV9Tuple(Y, X, h_, lat, lon, h, M, name, True)

    def reverse(self, enh_, n=None, h_=0, M=None, name=NN):
        # overloaded for the _ChLV.reverse.__doc__
        Y, X, h_, name = self._YXh_n4(enh_, n, h_, name=name)
        a, b, h = _ChLV._YXh_2abh3(Y, X, h_)
        a2, b2 = a**2, b**2

        lon = Fsum( 2.6779094, 4.728982 * a,
                               0.791484 * a * b,
                               0.1306   * a * b2,
                              -0.0436   * a * a2).fover(0.36)
        lat = Fsum(16.9023892, 3.238272 * b,
                              -0.270978 * a2,
                              -0.002528 * b2,
                              -0.0447   * a2 * b,
                              -0.014    * b2 * b).fover(0.36)
        return self._ChLV9Tuple(Y, X, h_, lat, lon, h, M, name, False)

    def _ChLV9Tuple(self, Y, X, h_, lat, lon, h, M, name, fw):
        '''(INTERNAL) Helper for C{.forward} and C{.reverse}.
        '''
        if M is not None:
            m =  self.forward if fw else self.reverse  # PYCHOK attr
            n = _DOT_(self.__class__.__name__, m.__name__)
            raise _NotImplementedError(unstr(n, M=M), txt=None)
        return ChLV9Tuple(Y, X, h_, lat, lon, h, self, self._t0, M, name=name)


class ChLVe(ChLVa):
    '''Conversion between I{WGS84 geodetic} and I{Swiss} projection coordinates
       using the U{Ellipsoidal approximate<https://www.SwissTopo.admin.CH/en/
       maps-data-online/calculation-services.html>} formulas, pp 10-11.

       @see: Older U{references<https://GitHub.com/alphasldiallo/Swisstopo-WGS84-LV03>}.
    '''
    __init__ = ChLVa.__init__
    '''New C{ChLVe} I{Ellipsoidal approximate WGS84-Swiss} converter, centered at I{Bern, Ch}.

       @see: L{ChLVa.__init__} for more information.
    '''

    def forward(self, latlonh, lon=None, height=0, M=None, name=NN):  # PYCHOK unused M
        # overloaded for the _ChLV.forward.__doc__
        lat, lon, h, name = _llhn4(latlonh, lon, height, name=name)
        a, b, h_ = _ChLV._llh2abh_3(lat, lon, h)
        F = Fhorner

        y1 = F(a,  0.2114285339, -0.010939608, -0.000002658, -0.00000853)
        y3 = F(a, -0.0000442327,  0.000004291, -0.000000309)
        y5 =       0.0000000197
        Y  = F(b,  0, y1, 0, y3, 0, y5).fover(_ChLV._YXm)  # + 2_600_000

        x0 = F(a,  0,             0.3087707463, 0.000075028, 0.000120435, 0, 0.00000007)
        x2 = F(a,  0.0037454089, -0.0001937927, 0.00000434, -0.000000376)
        x4 = F(a, -0.0000007346, 0.0000001444)
        X  = F(b, x0, 0, x2, 0, x4).fover(_ChLV._YXm)  # + 1_200_000

        return self._ChLV9Tuple(Y, X, h_, lat, lon, h, M, name, True)

    def reverse(self, enh_, n=None, h_=0, M=None, name=NN):  # PYCHOK unused M
        # overloaded for the _ChLV.reverse.__doc__
        Y, X, h_, name = self._YXh_n4(enh_, n, h_, name=name)
        a, b, h = _ChLV._YXh_2abh3(Y, X, h_)
        F = Fhorner

        a1  = F(b,  4.72973056, 0.7925714, 0.132812, 0.0255, 0.0048)
        a3  = F(b, -0.04427,   -0.0255,   -0.0096)
        a5  =       0.00096
        lon = F(a, 2.67825, a1, 0, a3, 0, a5).fsum()

        b0  = F(b, 16.902866,    3.23864877, -0.0025486, -0.013245, 0.000048)
        b2  = F(b, -0.27135379, -0.0450442,  -0.007553,  -0.00146)
        b4  = F(b,  0.002442,    0.00132)
        lat = F(a, b0, 0, b2, 0, b4).fsum()

        return self._ChLV9Tuple(Y, X, h_, lat, lon, h, M, name, False)


def tyr3d(tilt=INT0, yaw=INT0, roll=INT0, Vector=Vector3d, **Vector_kwds):
    '''Convert an attitude oriention into a (3-D) direction vector.

       @kwarg tilt: Pitch, elevation from horizontal (C{degrees}), negative down
                    (clockwise rotation along and around the x-axis).
       @kwarg yaw: Bearing, heading (compass C{degrees360}), clockwise from North
                   (counter-clockwise rotation along and around the z-axis).
       @kwarg roll: Roll, bank (C{degrees}), positive to the right and down
                    (clockwise rotation along and around the y-axis).

       @return: A named B{C{Vector}} instance or if B{C{Vector}} is C{None},
                a named L{Vector3Tuple}C{(x, y, z)}.

       @see: U{Yaw, pitch, and roll rotations<http://Planning.CS.UIUC.edu/node102.html>}
             and function L{pygeodesy.hartzell} argument C{los}.
    '''
    d = Attitude4Tuple(_0_0, tilt, yaw, roll).tyr3d
    return d if Vector is type(d) else (
           Vector3Tuple(d.x, d.y, d.z, name=d.name) if Vector is None else
                 Vector(d.x, d.y, d.z, **_xkwds(Vector_kwds, name=d.name)))  # PYCHOK indent


def _xLtp(ltp, *dflt):
    '''(INTERNAL) Validate B{C{ltp}}.
    '''
    if dflt and ltp is None:
        ltp = dflt[0]
    if isinstance(ltp, (LocalCartesian, Ltp)):
        return ltp
    raise _TypesError(_ltp_, ltp, Ltp, LocalCartesian)

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
