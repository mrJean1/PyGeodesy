
# -*- coding: utf-8 -*-

u'''Triaxal ellipsoid classes L{Triaxial} and I{unordered} L{Triaxial_} and Jacobi conformal projections
L{Conformal} and L{ConformalSphere}, transcoded from I{Karney}'s GeographicLib 2.5.2 C++ class U{JacobiConformal
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1JacobiConformal.html#details>} to pure
Python and miscellaneous classes L{BetaOmega2Tuple}, L{BetaOmega3Tuple} and L{Conformal2Tuple}, I{all kept
for backward copability}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2024) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib 2.5.2<https://GeographicLib.SourceForge.io>}
I{experimental} documentation.

@see: U{Geodesics on a triaxial ellipsoid<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
      Geodesics_on_a_triaxial_ellipsoid>} and U{Triaxial coordinate systems and their geometrical
      interpretation<https://OLD.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.

@var Triaxials.Amalthea: Triaxial(name='Amalthea', a=125000, b=73000, c=64000, e2ab=0.658944, e2bc=0.231375493, e2ac=0.737856, volume=2446253479595252, area=93239507787.490371704, area_p=93212299402.670425415)
@var Triaxials.Ariel: Triaxial(name='Ariel', a=581100, b=577900, c=577700, e2ab=0.01098327, e2bc=0.000692042, e2ac=0.011667711, volume=812633172614203904, area=4211301462766.580078125, area_p=4211301574065.829589844)
@var Triaxials.Earth: Triaxial(name='Earth', a=6378173.435, b=6378103.9, c=6356754.399999999, e2ab=0.000021804, e2bc=0.006683418, e2ac=0.006705077, volume=1083208241574987694080, area=510065911057441.0625, area_p=510065915922713.6875)
@var Triaxials.Enceladus: Triaxial(name='Enceladus', a=256600, b=251400, c=248300, e2ab=0.040119337, e2bc=0.024509841, e2ac=0.06364586, volume=67094551514082248, area=798618496278.596679688, area_p=798619018175.109985352)
@var Triaxials.Europa: Triaxial(name='Europa', a=1564130, b=1561230, c=1560930, e2ab=0.003704694, e2bc=0.000384275, e2ac=0.004087546, volume=15966575194402123776, area=30663773697323.51953125, area_p=30663773794562.45703125)
@var Triaxials.Io: Triaxial(name='Io', a=1829400, b=1819300, c=1815700, e2ab=0.011011391, e2bc=0.003953651, e2ac=0.014921506, volume=25313121117889765376, area=41691875849096.7421875, area_p=41691877397441.2109375)
@var Triaxials.Mars: Triaxial(name='Mars', a=3394600, b=3393300, c=3376300, e2ab=0.000765776, e2bc=0.009994646, e2ac=0.010752768, volume=162907283585817247744, area=144249140795107.4375, area_p=144249144150662.15625)
@var Triaxials.Mimas: Triaxial(name='Mimas', a=207400, b=196800, c=190600, e2ab=0.09960581, e2bc=0.062015624, e2ac=0.155444317, volume=32587072869017956, area=493855762247.691833496, area_p=493857714107.9375)
@var Triaxials.Miranda: Triaxial(name='Miranda', a=240400, b=234200, c=232900, e2ab=0.050915557, e2bc=0.011070811, e2ac=0.061422691, volume=54926187094835456, area=698880863325.757080078, area_p=698881306767.950317383)
@var Triaxials.Moon: Triaxial(name='Moon', a=1735550, b=1735324, c=1734898, e2ab=0.000260419, e2bc=0.000490914, e2ac=0.000751206, volume=21886698675223740416, area=37838824729886.09375, area_p=37838824733332.21875)
@var Triaxials.Tethys: Triaxial(name='Tethys', a=535600, b=528200, c=525800, e2ab=0.027441672, e2bc=0.009066821, e2ac=0.036259685, volume=623086233855821440, area=3528073490771.394042969, area_p=3528074261832.738769531)
@var Triaxials.WGS84_3: Triaxial(name='WGS84_3', a=6378171.36, b=6378101.609999999, c=6356751.84, e2ab=0.000021871, e2bc=0.006683505, e2ac=0.00670523, volume=1083207064030173855744, area=510065541435967.4375, area_p=510065546301413.5625)
@var Triaxials.WGS84_3r: Triaxial(name='WGS84_3r', a=6378172, b=6378102, c=6356752, e2ab=0.00002195, e2bc=0.006683577, e2ac=0.00670538, volume=1083207266220584468480, area=510065604942135.8125, area_p=510065609807745.0)
@var Triaxials.WGS84_35: Triaxial(name='WGS84_35', a=6378172, b=6378102, c=6356752.314245179, e2ab=0.00002195, e2bc=0.006683478, e2ac=0.006705281, volume=1083207319768789942272, area=510065621722018.125, area_p=510065626587483.3125)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.angles import _SinCos2,  Property_RO
from pygeodesy.basics import _isin, isLatLon
from pygeodesy.constants import EPS, EPS0, EPS02, _EPS2e4, INT0, \
                               _isfinite, isnear1, _over, _SQRT2_2, \
                               _0_0, _0_5, _1_0, _N_1_0, _64_0
from pygeodesy.datums import Datum, _spherical_datum, _WGS84,  _EWGS84, Fmt
# from pygeodesy.ellipsoids import Ellipsoid, _EWGS84  # from .datums
# from pygeodesy.elliptic import Elliptic  # _MODS
from pygeodesy.errors import _AssertionError, _ValueError, _xkwds_pop2
from pygeodesy.fmath import Fdot, fdot, hypot, hypot_,  fabs, sqrt
from pygeodesy.fsums import fsumf_, fsum1f_
from pygeodesy.interns import NN, _beta_, _distant_, _DMAIN_, _finite_, _height_, \
                             _inside_, _near_, _negative_, _not_, _null_, _opposite_, \
                             _outside_, _too_, _x_, _y_
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _lazyNamedEnumItem as _lazy, _name__, _NamedEnum, _Pass
from pygeodesy.namedTuples import LatLon3Tuple, _NamedTupleTo, Vector2Tuple, \
                                  Vector3Tuple, Vector4Tuple
# from pygeodesy.props import Property_RO  # from .triaxials.angles
# from pygeodesy.streprs import Fmt  # from .datums
from pygeodesy.triaxials.bases import Conformal5Tuple, _HeightINT0, _hypot2_1, \
                                     _not_ordered_, _OrderedTriaxialBase, _over0, \
                                     _otherV3d_, _over02, _sqrt0, TriaxialError, \
                                     _Triaxial3Base, _UnOrderedTriaxialBase
from pygeodesy.units import Degrees, Height_, Lat, Lon, Meter, Radians, Radius_, Scalar_
from pygeodesy.utily import atan2, atan2d, km2m, m2km
from pygeodesy.vector3d import _otherV3d, Vector3d

# from math import fabs, sqrt  # from .fmath

__all__ = _ALL_LAZY.triaxials_triaxial5
__version__ = '25.11.29'

_omega_ = 'omega'
_TRIPS  =  359  # Eberly 1074?


class _NamedTupleToX(_NamedTupleTo):  # in .testNamedTuples
    '''(INTERNAL) Base class for L{BetaOmega2Tuple},
       L{BetaOmega3Tuple} and L{Conformal2Tuple}.
    '''
    def _toDegrees(self, name, **toDMS_kwds):
        '''(INTERNAL) Convert C{self[0:2]} to L{Degrees} or C{toDMS}.
        '''
        return self._toX3U(_NamedTupleTo._Degrees3, Degrees, name, *self, **toDMS_kwds)

    def _toRadians(self, name):
        '''(INTERNAL) Convert C{self[0:2]} to L{Radians}.
        '''
        return self._toX3U(_NamedTupleTo._Radians3, Radians, name, *self)

    def _toX3U(self, _X3, U, name, a, b, *c, **kwds):
        a, b, s = _X3(self, a, b, **kwds)
        if s is None or name:
            n = self._name__(name)
            s = self.classof(a, b, *c, name=n).reUnit(U, U).toUnits()
        return s


class BetaOmega2Tuple(_NamedTupleToX):
    '''2-Tuple C{(beta, omega)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in L{Radians} (or
       L{Degrees}).
    '''
    _Names_ = (_beta_, _omega_)
    _Units_ = (_Pass,  _Pass)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{BetaOmega2Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with C{beta} and
                    C{omega} both in L{Degrees} or as L{toDMS} strings
                    provided some B{C{toDMS_kwds}} keyword arguments are
                    specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{BetaOmega2Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with C{beta} and C{omega}
                    both in L{Radians}.
        '''
        return self._toRadians(name)


class BetaOmega3Tuple(_NamedTupleToX):
    '''3-Tuple C{(beta, omega, height)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in L{Radians} (or L{Degrees})
       and the C{height}, rather the (signed) I{distance} to the triaxial's
       surface (measured along the radial line to the triaxial's center)
       in C{meter}, conventionally.
    '''
    _Names_ = BetaOmega2Tuple._Names_ + (_height_,)
    _Units_ = BetaOmega2Tuple._Units_ + ( Meter,)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{BetaOmega3Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with
                    C{beta} and C{omega} both in L{Degrees} or as
                    L{toDMS} strings provided some B{C{toDMS_kwds}}
                    keyword arguments are specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{BetaOmega3Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding this name.

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta}
                    and C{omega} both in L{Radians}.
        '''
        return self._toRadians(name)

    def to2Tuple(self, **name):
        '''Reduce this L{BetaOmega3Tuple} to a L{BetaOmega2Tuple}.

           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding this name.
        '''
        return BetaOmega2Tuple(*self[:2], name=self._name__(name))


class Conformal2Tuple(_NamedTupleToX):
    '''2-Tuple C{(x, y)} with a I{Jacobi Conformal} C{x} and C{y}
       projection, both in L{Radians} (or L{Degrees}).
    '''
    _Names_ = (_x_,   _y_)
    _Units_ = (_Pass, _Pass)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{Conformal2Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Conformal2Tuple}C{(x, y)} with C{x} and C{y} both
                    in L{Degrees} or as L{toDMS} strings provided some
                    B{C{toDMS_kwds}} keyword arguments are specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{Conformal2Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Conformal2Tuple}C{(x, y)} with C{x} and C{y} both in L{Radians}.
        '''
        return self._toRadians(name)

    def to5Tuple(self, b_conformal, **z_scale_name):
        '''Return this L{Conformal2Tuple} as a L{Conformal5Tuple}.

           @arg b_conformal: Middle semi-axis (C{meter}, conventionally)
                             or the original L{Conformal} of this 2-tuple.
           @kwarg z_scale_name: Optional C{B{z}=0} meter, C{B{scale}=NAN}
                                and C{B{name}=NN} (C{str}).

           @return: A L{Conformal5Tuple}.
        '''
        if isinstance(b_conformal, Conformal):
            b = b_conformal.b
        else:
            b = Radius_(b=b_conformal)
        x, y =  self.toRadians()
        x, y = _over(x, b), _over(y, b)
        return Conformal5Tuple(x, y, **z_scale_name)


class Triaxial_(_UnOrderedTriaxialBase):
    '''I{Unordered} triaxial ellipsoid.

       Triaxial ellipsoids with right-handed semi-axes C{a}, C{b} and C{c}, oriented
       such that the large principal ellipse C{ab} is the equator I{Z}=0, I{beta}=0,
       while the small principal ellipse C{ac} is the prime meridian, plane I{Y}=0,
       I{omega}=0.

       The four umbilic points, C{abs}(I{omega}) = C{abs}(I{beta}) = C{PI/2}, lie on
       the middle principal ellipse C{bc} in plane I{X}=0, I{omega}=C{PI/2}.

       @note: I{Geodetic} C{lat}- and C{lon}gitudes are in C{degrees}, I{geodetic}
              C{phi} and C{lam}bda are in C{radians}, but I{ellipsoidal} lat- and
              longitude C{beta} and C{omega} are in L{Radians} by default (or in
              L{Degrees} if converted).
    '''
    if _FOR_DOCS:
        __init__ = _UnOrderedTriaxialBase.__init__


class Triaxial(_OrderedTriaxialBase):
    '''I{Ordered} triaxial ellipsoid.

       @see: L{Triaxial_} for more information.
    '''
    if _FOR_DOCS:
        __init__ = _OrderedTriaxialBase.__init__

    def forwardBetaOmega(self, beta, omega, height=0, **unit_name):
        '''Convert I{ellipsoidal} lat- C{beta}, longitude C{omega} and C{height}
           to cartesian.

           @arg beta: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omega: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg height: Height above or below the triaxial's surface (C{meter},
                          same units as this triaxial's semi-axes.
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar
                             C{B{unit}=}L{Radians} (or L{Degrees}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseBetaOmega} and U{equations (23-25)<https://
                 OLD.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        unit, name = _xkwds_pop2(unit_name, unit=Radians)
        if height:
            z = self._Height(height) + self.c
            if z > 0:
                z2 = z**2
                x  = z * _sqrt0(_1_0 + self._a2c2 / z2)
                y  = z * _sqrt0(_1_0 + self._b2c2 / z2)
            else:
                x  = y = z = _0_0
        else:
            x, y, z = self._abc3
        if z:  # and x and y:
            sa, ca = _SinCos2(beta,  unit)
            sb, cb = _SinCos2(omega, unit)

            r  = self._a2b2_a2c2
            x *= cb * (_sqrt0(ca**2 + sa**2 * r) if r else fabs(ca))
            y *= ca *   sb
            z *= sa * (_sqrt0(_1_0  - cb**2 * r) if r else _1_0)
        return Vector3Tuple(x, y, z, **name)

    def forwardBetaOmega_(self, sbeta, cbeta, somega, comega, **name):
        '''Convert I{ellipsoidal} lat- and longitude C{beta} and C{omega}
           to cartesian coordinates I{on this triaxial's surface}.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)} on the surface.

           @raise TriaxialError: This triaxial is near-spherical.

           @see: Method L{Triaxial.reverseBetaOmega}, U{Triaxial ellipsoid coordinate
                 system<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
                 Triaxial_ellipsoid_coordinate_system>} and U{equations (23-25)<https://
                 OLD.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        t = self._radialTo3(sbeta, cbeta, somega, comega)
        return Vector3Tuple(*t, **name)

    def forwardCartesian(self, x_xyz, y=None, z=None, normal=True, **eps_name):
        '''Project any cartesian to a cartesian I{on this triaxial's surface}.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg normal: If C{True}, the projection is C{perpendicular} to the surface,
                          otherwise C{radial} to the center of this triaxial (C{bool}).
           @kwarg eps_name: Root finder tolerance C{B{eps}=EPS} and optional
                            C{B{name}="height4"} (C{str}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)}.

           @see: Method L{Triaxial.reverseCartesian} to reverse the projection and
                 function L{height4<triaxials.triaxial5.height4>} for more details.
        '''
        return self.height4(x_xyz, y, z, normal=normal, **eps_name)

    def forwardLatLon(self, lat, lon, height=0, **unit_name):
        '''Convert I{geodetic} lat-, longitude and height to cartesian.

           @arg lat: Geodetic latitude (C{Ang} or B{C{unit}}).
           @arg lon: Geodetic longitude (C{Ang} or B{C{unit}}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's semi-axes).
           @kwarg unit_name: Optional scalar C{B{unit}=}L{Degrees} and
                       C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseLatLon} and U{equations (9-11)<https://
                 OLD.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        unit, name = _xkwds_pop2(unit_name, unit=Degrees)
        return self._forwardLatLon3(height, name, *(_SinCos2(lat, unit) +
                                                    _SinCos2(lon, unit)))

    def forwardLatLon_(self, slat, clat, slon, clon, height=0, **name):
        '''Convert I{geodetic} lat-, longitude and height to cartesian.

           @arg slat: Geodetic latitude C{sin(lat)} (C{scalar}).
           @arg clat: Geodetic latitude C{cos(lat)} (C{scalar}).
           @arg slon: Geodetic longitude C{sin(lon)} (C{scalar}).
           @arg clon: Geodetic longitude C{cos(lon)} (C{scalar}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's semi-axes).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseLatLon} and U{equations (9-11)<https://
                 OLD.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca = self._norm2(slat, clat)
        sb, cb = self._norm2(slon, clon)
        return self._forwardLatLon3(height, name, sa, ca, sb, cb)

    def _forwardLatLon3(self, height, name, sa, ca, sb, cb):  # name always **name
        '''(INTERNAL) Helper for C{.forwardLatLon} and C{.forwardLatLon_}.
        '''
        h = self._Height(height)
        x = ca * cb
        y = ca * sb
        z = sa
        # 1 - (1 - (c/a)**2) * sa**2 - (1 - (b/a)**2) * ca**2 * sb**2
        t = fsumf_(_1_0, -self.e2ac * z**2, -self.e2ab * y**2)
        n = self.a / _sqrt0(t)  # prime vertical
        x *= h + n
        y *= h + n * self._b2_a2
        z *= h + n * self._c2_a2
        return Vector3Tuple(x, y, z, **name)

    def _Height(self, height):
        '''(INTERNAL) Validate a C{height}.
        '''
        return Height_(height=height, low=-self.c, Error=TriaxialError)

    def reverseBetaOmega(self, x_xyz, y=None, z=None, **name):
        '''Convert cartesian to I{ellipsoidal} lat- and longitude, C{beta}, C{omega}
           and height.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta} and
                    C{omega} in L{Radians} and (radial) C{height} in C{meter}, same
                    units as this triaxial's semi-axes.

           @see: Methods L{Triaxial.forwardBetaOmega} and L{Triaxial.forwardBetaOmega_}
                 and U{equations (21-22)<https://OLD.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        a, b, h = _reverseLatLon3(v, atan2, v, self.forwardBetaOmega_)
        return BetaOmega3Tuple(Radians(beta=a), Radians(omega=b), h, **name)

    def reverseCartesian(self, x_xyz, y=None, z=None, height=0, normal=True, eps=_EPS2e4, **name):
        '''"Unproject" a cartesian I{off} this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                       ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg height: Height above or below this triaxial's surface (C{meter}, same
                          units as this triaxial's semi-axes).
           @kwarg normal: If C{True}, B{C{height}} is C{perpendicular} to the surface,
                          otherwise C{radial} to the center of this triaxial (C{bool}).
           @kwarg eps: Tolerance for on-surface test (C{scalar}), see method
                       L{Triaxial.sideOf}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @raise TrialError: Cartesian B{C{x_xyz}} or C{(x, y, z)} not on this triaxial's
                              surface.

           @see: Methods L{Triaxial.forwardCartesian} and L{Triaxial.height4}.
        '''
        h, name = _xkwds_pop2(name, h=height)  # h=height for backward compatibility
        v = _otherV3d_(x_xyz, y, z, **name)
        _ =  self._sideOn(v, eps=eps)
        h = _HeightINT0(h)
        if h:
            if normal:
                v = v.plus(self.normal3d(*v.xyz, length=h))
            elif v.length > EPS0:
                v = v.times(_1_0 + (h / v.length))
        return v.xyz  # Vector3Tuple

    def reverseLatLon(self, x_xyz, y=None, z=None, **name):
        '''Convert cartesian to I{geodetic} lat-, longitude and height.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)} with C{lat} and C{lon}
                    in C{degrees} and (radial) C{height} in C{meter}, same units
                    as this triaxial's semi-axes.

           @see: Methods L{Triaxial.forwardLatLon} and L{Triaxial.forwardLatLon_}
                 and U{equations (4-5)<https://OLD.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        s =  v.times_(self._c2_a2,  # == 1 - e_sub_x**2
                      self._c2_b2,  # == 1 - e_sub_y**2
                     _1_0)
        a, b, h = _reverseLatLon3(s, atan2d, v, self.forwardLatLon_)
        return LatLon3Tuple(Lat(a), Lon(b), h, **name)


class TriaxialB(_Triaxial3Base):
    '''I{Ordered} triaxial ellipsoid specified by its middle semi-axis and shape.

       @see: L{Triaxial} for details and more information.
    '''
    def __init__(self, b, e2=_0_0, k2=_1_0, kp2=_0_0, **name):
        '''New L{TriaxialB} triaxial.

           @arg b: Middle semi-axis (C{meter}, conventionally).
           @kwarg e2: Excentricty I{squared} (C{scalar, e2 = (a**2 - c**2) / B{b}**2}).
           @kwarg k2: Oblateness I{squared} (C{scalar, k2 = (C{b}**2 - c**2) / (a**2 - c**2)}).
           @kwarg kp2: Prolateness I{squared} (C{scalar, kp2 = (a**2 - C{b}**2) / (a**2 - c**2)}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @note: The semi-axes must be ordered as C{B{a} >= B{b} >= B{c} > 0} but
                  may be spherical, C{B{e2} == 0}, i.e. C{B{a} == B{c}}.  In that case
                  C{B{k2} == 0} represents a I{prolate sphere}, C{B{k2} == 1} an
                  I{oblate sphere}, otherwise a I{triaxial sphere} with C{0 < B{k2} < 1}.

           @note: The B{C{k2}} and B{C{kp2}} arguments are normalized, C{B{k2} + B{kp2} == 1}.

           @raise TriaxialError: Semi-axes unordered or invalid.
        '''
        self._init_abc3_e2_k2_kp2(Radius_(b=b), e2, k2, kp2, **name)


class Conformal(Triaxial):
    '''This is a I{Jacobi Conformal} projection of a triaxial ellipsoid to a plane where
       the C{X} and C{Y} grid lines are straight.

       I{Ellipsoidal} coordinates I{beta} and I{omega} are converted to Jacobi Conformal
       I{y} respectively I{x} separately.  Jacobi's coordinates have been multiplied
       by C{sqrt(B{a}**2 - B{c}**2) / (2 * B{b})} so that the customary results are
       returned in the case of an ellipsoid of revolution.

       Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2014-2024) and
       licensed under the MIT/X11 License.

       @note: This constructor can I{not be used to specify a sphere}, see alternate
              L{ConformalSphere}.

       @see: L{Triaxial}, C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1JacobiConformal.html#details>}, U{Jacobi's conformal
             projection<https://GeographicLib.SourceForge.io/C++/doc/jacobi.html>} and Jacobi,
             C. G. J. I{U{Vorlesungen über Dynamik<https://Books.Google.com/books?
             id=ryEOAAAAQAAJ&pg=PA212>}}, page 212ff.
    '''
    if _FOR_DOCS:
        __init__ = Triaxial.__init__

    @Property_RO
    def _a2_b(self):
        return self._a2_b2 * self.b

    @Property_RO
    def _c2_b(self):
        return self._c2_b2 * self.b

    def x(self, omega, unit=Radians):
        '''Compute a I{Jacobi Conformal} C{x} projection.

           @arg omega: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg unit: Unit of scalar B{C{omega}} (or C{Degrees}).

           @return: The C{x} projection (L{Meter}), same units as
                    this triaxial's semi-axes.
        '''
        s, c = _SinCos2(omega, unit)
        return Meter(x=self._x(s, c, self._a2_b))

    def _x(self, s, c, a2_b_):
        '''(INTERNAL) Helper for C{.x}, C{.xR_} and C{.xy}.
        '''
        s, c = self._norm2(s, c, self.a)
        return self._xE.fPi(s, c) * a2_b_

    @Property_RO
    def _xE(self):
        '''(INTERNAL) Get the x-elliptic function.
        '''
        k2, kp2 = self._k2E_kp2E
        # -a2b2 / b2 == (b2 - a2) / b2 == 1 - a2 / b2 == 1 - a2_b2
        return self._Elliptic(k2, _1_0 - self._a2_b2, kp2, self._a2_b2)

    def xR(self, omega, unit=Radians):
        '''Compute a I{Jacobi Conformal} C{x} projection.

           @arg omega: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg unit: Unit of scalar B{C{omega}} (or C{Degrees}).

           @return: The C{x} projection (L{Radians}).
        '''
        return self.xR_(*_SinCos2(omega, unit))

    def xR_(self, somega, comega):
        '''Compute a I{Jacobi Conformal} C{x} projection.

           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).

           @return: The C{x} projection (L{Radians}).
        '''
        return Radians(x=self._x(somega, comega, self._a2_b2))

    def xy(self, beta, omega, **unit_name):
        '''Compute a I{Jacobi Conformal} C{x} and C{y} projection.

           @arg beta: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omega: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg unit_name: Optional scalar C{B{unit}=}L{Radians} and
                       name (C{str}), overriding C{B{name}="xyR2"}.

           @return: A (L{Vector2Tuple}C{(x, y)}), both in C{Meter},
                    same units as this triaxial's semi-axes..
        '''
        unit, name = _xkwds_pop2(unit_name, unit=Radians)
        return Vector2Tuple(self.x(omega, unit=unit),
                            self.y(beta,  unit=unit),
                            name=_name__(name, name__=self.xy))

    @Property_RO
    def xyQ2(self):
        '''Get the I{Jacobi Conformal} quadrant size in C{meter}
           (L{Vector2Tuple}C{(x, y)}).
        '''
        return Vector2Tuple(Meter(x=self._a2_b * self._xE.cPi),
                            Meter(y=self._c2_b * self._yE.cPi),
                            name=Conformal.xyQ2.name)

    @Property_RO
    def xyQR2(self):
        '''Get the I{Jacobi Conformal} quadrant size in C{Radians}
           (L{Conformal2Tuple}C{(x, y)}).
        '''
        return Conformal2Tuple(Radians(x=self._a2_b2 * self._xE.cPi),
                               Radians(y=self._c2_b2 * self._yE.cPi),
                               name=Conformal.xyQR2.name)

    def xyR2(self, beta, omega, **unit_name):
        '''Compute a I{Jacobi Conformal} C{x} and C{y} projection.

           @arg beta: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omega: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg unit_name: Optional scalar C{B{unit}=}L{Radians} and
                       name (C{str}), overriding C{B{name}="xyR2"}.

           @return: A L{Conformal2Tuple}C{(x, y)}, both in C{Radians}.
        '''
        unit, name  = _xkwds_pop2(unit_name, unit=Radians)
        sb_cb_so_co = _SinCos2(beta, unit) + _SinCos2(omega, unit)
        return self.xyR2_(*sb_cb_so_co, name=_name__(name, name__=self.xyR2))

    def xyR2_(self, sbeta, cbeta, somega, comega, **name):
        '''Compute a I{Jacobi Conformal} C{x} and C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional name (C{str}), overriding C{B{name}="xyR2_"}.

           @return: A L{Conformal2Tuple}C{(x, y)}.
        '''
        return Conformal2Tuple(self.xR_(somega, comega),
                               self.yR_(sbeta,  cbeta),
                               name=_name__(name, name__=self.xyR2_))

    def y(self, beta, unit=Radians):
        '''Compute a I{Jacobi Conformal} C{y} projection.

           @arg beta: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @kwarg unit: Unit of scalar B{C{beta}} (or C{Degrees}).

           @return: The C{y} projection (L{Meter}), same units as
                    this triaxial's semi-axes.
        '''
        s, c = _SinCos2(beta, unit)
        return Meter(y=self._y(s, c, self._c2_b))

    def _y(self, s, c, c2_b_):
        '''(INTERNAL) Helper for C{.y}, C{.yR_} and C{.xy}.
        '''
        s, c = self._norm2(s, c, self.c)
        return self._yE.fPi(s, c) * c2_b_

    @Property_RO
    def _yE(self):
        '''(INTERNAL) Get the y-elliptic function.
        '''
        k2, kp2 = self._k2E_kp2E
        # b2c2 / b2 == (b2 - c2) / b2 == 1 - c2 / b2 == e2bc
        return self._Elliptic(kp2, self.e2bc, k2, self._c2_b2)

    def yR(self, beta, unit=Radians):
        '''Compute a I{Jacobi Conformal} C{y} projection.

           @arg beta: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @kwarg unit: Unit of scalar B{C{beta}} (or C{Degrees}).

           @return: The C{y} projection (L{Radians}).
        '''
        return self.yR_(*_SinCos2(beta, unit))

    def yR_(self, sbeta, cbeta):
        '''Compute a I{Jacobi Conformal} C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).

           @return: The C{y} projection (L{Radians}).
        '''
        return Radians(y=self._y(sbeta, cbeta, self._c2_b2))


class ConformalSphere(Conformal):
    '''Alternate, I{Jacobi Conformal projection} on a I{spherical} triaxial.

       @see: L{Conformal<triaxials.triaxial5.Conformal>} for more information.
    '''
    _ab = _bc = 0

    def __init__(self, radius_conformal, ab=0, bc=0, **name):
        '''New L{ConformalSphere}.

           @arg radius_conformal: Radius (C{scalar}, conventionally in C{meter})
                       or an other L{ConformalSphere} or L{Conformal}.
           @kwarg ab: Relative magnitude of C{B{a} - B{b}} (C{meter}, same units
                      as C{scalar B{radius}}.
           @kwarg bc: Relative magnitude of C{B{b} - B{c}} (C{meter}, same units
                      as C{scalar B{radius}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TriaxialError: Invalid B{C{radius_conformal}}, negative B{C{ab}},
                                 negative B{C{bc}} or C{(B{ab} + B{bc})} not positive.

           @note: If B{C{radius_conformal}} is a L{ConformalSphere} and if B{C{ab}}
                  and B{C{bc}} are both zero or C{None}, the B{C{radius_conformal}}'s
                  C{ab}, C{bc}, C{a}, C{b} and C{c} are copied.
        '''
        try:
            r = radius_conformal
            if isinstance(r, Triaxial):  # ordered only
                t = r._abc3
                j = isinstance(r, ConformalSphere) and not bool(ab or bc)
            else:
                t = (Radius_(radius=r),) * 3
                j =  False
            self._ab = r.ab if j else Scalar_(ab=ab)  # low=0
            self._bc = r.bc if j else Scalar_(bc=bc)  # low=0
            if (self.ab + self.bc) <= 0:
                raise ValueError(_negative_)
            a, _, c = self._abc3 = t
            if not (a >= c and _isfinite(self._a2b2)
                           and _isfinite(self._a2c2)):
                raise ValueError(_not_(_finite_))
        except (TypeError, ValueError) as x:
            raise TriaxialError(radius=r, ab=ab, bc=bc, cause=x)
        if name:
            self.name = name

    @Property_RO
    def ab(self):
        '''Get relative magnitude C{a - b} (C{meter}, same units as B{C{a}}).
        '''
        return self._ab

    @Property_RO
    def _a2b2(self):
        '''(INTERNAL) Get C{a**2 - b**2} == ab * (a + b).
        '''
        a, b, _ = self._abc3
        return self.ab * (a + b)

    @Property_RO
    def _a2c2(self):
        '''(INTERNAL) Get C{a**2 - c**2} == a2b2 + b2c2.
        '''
        return self._a2b2 + self._b2c2

    @Property_RO
    def bc(self):
        '''Get relative magnitude C{b - c} (C{meter}, same units as B{C{a}}).
        '''
        return self._bc

    @Property_RO
    def _b2c2(self):
        '''(INTERNAL) Get C{b**2 - c**2} == bc * (b + c).
        '''
        _, b, c = self._abc3
        return self.bc * (b + c)

    @Property_RO
    def radius(self):
        '''Get radius (C{meter}, conventionally).
        '''
        return self.a


class Triaxials(_NamedEnum):
    '''(INTERNAL) L{Triaxial} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, *abc, **name):
        '''(INTERNAL) Instantiate the C{Triaxial}.
        '''
        a, b, c = map(km2m, abc)
        return Triaxial(a, b, c, **name)

Triaxials = Triaxials(Triaxial, Triaxial_)  # PYCHOK singleton
'''Some pre-defined L{Triaxial}s, all I{lazily} instantiated.'''
# <https://ArxIV.org/pdf/1909.06452.pdf> Table 1 Semi-axes in Km
# <https://www.JPS.NASA.gov/education/images/pdf/ss-moons.pdf>
# <https://link.Springer.com/article/10.1007/s00190-022-01650-9>
# <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Constants.html>
_abc84_35 = (_EWGS84.a + 35), (_EWGS84.a - 35), _EWGS84.b
Triaxials._assert(                 # a (Km)       b (Km)      c (Km)       planet
    Amalthea  = _lazy('Amalthea',  125.0,        73.0,      _64_0),      # Jupiter
    Ariel     = _lazy('Ariel',     581.1,       577.9,      577.7),      # Uranus
    Earth     = _lazy('Earth',    6378.173435, 6378.1039,  6356.7544),
    Enceladus = _lazy('Enceladus', 256.6,       251.4,      248.3),      # Saturn
    Europa    = _lazy('Europa',   1564.13,     1561.23,    1560.93),     # Jupiter
    Io        = _lazy('Io',       1829.4,      1819.3,     1815.7),      # Jupiter
    Mars      = _lazy('Mars',     3394.6,      3393.3,     3376.3),
    Mimas     = _lazy('Mimas',     207.4,       196.8,      190.6),      # Saturn
    Miranda   = _lazy('Miranda',   240.4,       234.2,      232.9),      # Uranus
    Moon      = _lazy('Moon',     1735.55,     1735.324,   1734.898),    # Earth
    Tethys    = _lazy('Tethys',    535.6,       528.2,      525.8),      # Saturn
    WGS84_3   = _lazy('WGS84_3',  6378.17136,  6378.10161, 6356.75184),  # C++
    WGS84_3r  = _lazy('WGS84_3r', 6378.172,    6378.102,   6356.752),    # C++, rounded
    WGS84_35  = _lazy('WGS84_35', *map(m2km, _abc84_35)))
del _abc84_35, _EWGS84


def _hartzell3(pov, los, Tun):  # in .Ellipsoid.hartzell4, .formy.hartzell
    '''(INTERNAL) Hartzell's "Satellite Line-of-Sight Intersection ...",
       formula from a Point-Of-View to an I{un-/ordered} Triaxial.
    '''
    def _toUvwV3d(los, pov):
        try:  # pov must be LatLon or Cartesian if los is a Los
            v = los.toUvw(pov)
        except (AttributeError, TypeError):
            v = _otherV3d(los=los)
        return v

    p3 = _otherV3d(pov=pov.toCartesian() if isLatLon(pov) else pov)
    if los is True:  # normal
        a, b, c, d, i = _plumbTo5(p3.x, p3.y, p3.z, Tun)
        return type(p3)(a, b, c), d, i

    u3 = p3.negate() if los is False or los is None else _toUvwV3d(los, pov)

    a, b, c, T = Tun._ordered4

    a2     =  T.a2  # largest, factored out
    b2, p2 = (T.b2, T._b2_a2) if b != a else (a2, _1_0)
    c2, q2 = (T.c2, T._c2_a2) if c != a else (a2, _1_0)

    p3 = T._order3d(p3)
    u3 = T._order3d(u3).unit()  # unit vector, opposing signs

    x2, y2, z2 = p3.x2y2z23  # p3.times_(p3).xyz3
    ux, vy, wz = u3.times_(p3).xyz3
    u2, v2, w2 = u3.x2y2z23  # u3.times_(u3).xyz3

    t = (p2 * c2),  c2, b2
    m = fdot(t, u2, v2, w2)  # a2 factored out
    if m < EPS0:  # zero or near-null LOS vector
        raise _ValueError(_near_(_null_))

    r = fsumf_(b2 * w2,      c2 * v2,      -v2 * z2,      vy * wz * 2,
              -w2 * y2,     -u2 * y2 * q2, -u2 * z2 * p2, ux * wz * 2 * p2,
              -w2 * x2 * p2, b2 * u2 * q2, -v2 * x2 * q2, ux * vy * 2 * q2)
    if r > 0:  # a2 factored out
        r = sqrt(r) * b * c  # == a * a * b * c / a2
    elif r < 0:  # LOS pointing away from or missing the triaxial
        raise _ValueError(_opposite_ if max(ux, vy, wz) > 0 else _outside_)

    d = Fdot(t, ux, vy, wz).fadd_(r).fover(m)  # -r for antipode, a2 factored out
    if d > 0:  # POV inside or LOS outside or missing the triaxial
        s = fsumf_(_N_1_0, _over(x2, a2), _over(y2, b2), _over(z2, c2))  # like _sideOf
        raise _ValueError(_outside_ if s > 0 else _inside_)
    elif fsum1f_(x2, y2, z2, -d*d) < 0:  # d past triaxial's center
        raise _ValueError(_too_(_distant_))

    v = p3.minus(u3.times(d))  # cartesian type(pov) or Vector3d
    h = p3.minus(v).length  # distance to pov == -d
    return T._order3d(v, reverse=True), h, None


def hartzell4(pov, los=False, tri_biax=_WGS84, **name):
    '''Compute the intersection of a tri-/biaxial ellipsoid and a Line-Of-Sight from
       a Point-Of-View outside.

       @arg pov: Point-Of-View outside the tri-/biaxial (C{Cartesian}, L{Ecef9Tuple},
                 C{LatLon} or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to the tri-/biaxial (L{Los}, L{Vector3d})
                   or C{True} for the I{normal, perpendicular, plumb} to the surface of
                   the tri-/biaxial or C{False} or C{None} to point to its center.
       @kwarg tri_biax: A triaxial (L{Triaxial}, L{Triaxial_}, L{Conformal} or
                        L{ConformalSphere}) or biaxial ellipsoid (L{Datum},
                        L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple}) or spherical earth
                        radius (C{scalar}, conventionally in C{meter}).
       @kwarg name: Optional name (C{str}), overriding C{B{name}="hartzell4"}.

       @return: L{Vector4Tuple}C{(x, y, z, h)} on the tri-/biaxial's surface, with C{h}
                the distance from B{C{pov}} to C{(x, y, z)} I{along the} B{C{los}}, all
                in C{meter}, conventionally.

       @raise TriaxialError: Invalid B{C{pov}} or B{C{pov}} inside the tri-/biaxial or
                             invalid B{C{los}} or B{C{los}} points outside or away from
                             the tri-/biaxial.

       @raise TypeError: Invalid B{C{tri_biax}}, C{ellipsoid} or C{datum}.

       @see: Class L{pygeodesy3.Los}, functions L{pygeodesy.tyr3d} and L{pygeodesy.hartzell}
             and U{lookAtSpheroid<https://PyPI.org/project/pymap3d>} and U{"Satellite
             Line-of-Sight Intersection with Earth"<https://StephenHartzell.Medium.com/
             satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}.
    '''
    T = _tri_biaxial(tri_biax, hartzell4)
    try:
        v, h, i = _hartzell3(pov, los, T)
    except Exception as x:
        raise TriaxialError(pov=pov, los=los, tri_biax=tri_biax, cause=x)
    n = _name__(name, name__=hartzell4)  # typename
    return Vector4Tuple(v.x, v.y, v.z, h, iteration=i, name=n)


def height4(x_xyz, y=None, z=None, tri_biax=_WGS84, normal=True, eps=EPS, **name):
    '''Compute the projection on and the height above or below a tri-/biaxial ellipsoid's surface.

       @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian}, L{Ecef9Tuple},
                   L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
       @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}, ignored
                 otherwise.
       @kwarg z: Z component (C{scalar}), like B{C{y}}.
       @kwarg normal: If C{True}, the projection is the I{perpendicular, plumb} to the
                      tri-/biaxial's surface, otherwise the C{radial} line to the center
                      of the tri-/biaxial (C{bool}).
       @kwarg eps: Tolerance for root finding and validation (C{scalar}), use a negative
                   value to skip validation.
       @kwarg name: Optional C{B{name}="height4"} (C{str}).

       @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates C{x}, C{y}
                and C{z} of the projection on or the intersection with and with the height
                C{h} above or below the tri-/biaxial's surface in C{meter}, conventionally.

       @raise TriaxialError: Non-cartesian B{C{xyz}}, invalid B{C{eps}}, no convergence in
              root finding or validation failed.

       @see: Methods L{Triaxial.normal3d} and L{Ellipsoid.height4}, I{Eberly}'s U{Distance from a Point to ...
             <https://www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>} and I{Bektas}'
             U{Shortest Distance from a Point to Triaxial Ellipsoid<https://www.ResearchGate.net/publication/
             272149005_SHORTEST_DISTANCE_FROM_A_POINT_TO_TRIAXIAL_ELLIPSOID>}.
    '''
    v = _otherV3d_(x_xyz, y, z)
    T = _tri_biaxial(tri_biax, height4)
    r =  T.isSpherical

    i, h = None, v.length
    if h < EPS0:  # EPS
        x = y = z = _0_0
        h -= min(T._abc3)  # nearest
    elif r:  # .isSpherical
        x, y, z = v.times(r / h).xyz3
        h -= r
    else:
        x, y, z = v.xyz3
        try:
            if normal:  # plumb to surface
                x, y, z, h, i = _plumbTo5(x, y, z, T, eps=eps)
            else:  # radial to center
                x, y, z = T._radialTo3(z, hypot(x, y), y, x)
                h = v.minus_(x, y, z).length
        except Exception as e:
            raise TriaxialError(x=x, y=y, z=z, cause=e)
        if h > 0 and T.sideOf(v, eps=EPS0) < 0:
            h = -h  # inside
    n = _name__(name, name__=height4)  # typename
    return Vector4Tuple(x, y, z, h, iteration=i, name=n)


def _ordered(a, b, c):
    '''(INTERNAL) Is C{a >= b >= c > 0}?
    '''
    if not (_isfinite(a) and a >= b >= c > 0):
        raise TriaxialError(a=a, b=b, c=c, txt=_not_ordered_)


def _plumbTo3(px, py, E, eps=EPS):  # in .ellipsoids.Ellipsoid.height4
    '''(INTERNAL) Nearest point in 1st quadrant on a 2-D ellipse.
    '''
    a, b = E.a, E.b
    if min(px, py, a, b) < EPS0:
        raise _AssertionError(px=px, py=py, a=a, b=b, E=E)

    a2 = a - b * E.b_a
    b2 = b - a * E.a_b
    tx = ty = _SQRT2_2
    for i in range(16):  # max 5
        ex = tx**3 * a2
        ey = ty**3 * b2

        qx = px - ex
        qy = py - ey
        q  = hypot(qx, qy)
        if q < EPS0:  # PYCHOK no cover
            break
        r = hypot(ex - tx * a,
                  ey - ty * b) / q

        sx, tx = tx, min(_1_0, max(0, (ex + qx * r) / a))
        sy, ty = ty, min(_1_0, max(0, (ey + qy * r) / b))
        t = hypot(ty, tx)
        if t < EPS0:  # PYCHOK no cover
            break
        tx = tx / t  # /= chokes PyChecker
        ty = ty / t
        if fabs(sx - tx) < eps and \
           fabs(sy - ty) < eps:
            break

    tx *= a / px
    ty *= b / py
    return tx, ty, i  # x and y as fractions


def _plumbTo4(x, y, a, b, eps=EPS):
    '''(INTERNAL) Nearest point on and distance to a 2-D ellipse, I{unordered}.

       @see: Function C{_plumbTo3} and I{Eberly}'s U{Distance from a Point to ... an Ellipse ...
             <https://www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    if b > a:
        b, a, d, i = _plumbTo4(y, x, b, a, eps=eps)
        return a, b, d, i

    if not (b > 0 and _isfinite(a)):
        raise _ValueError(a=a, b=b)

    i = None
    if y:
        if x:
            u =  fabs(x / a)
            w =  fabs(y / b)
            g = _hypot2_1(u, w)
            if fabs(g) > EPS02:
                r = (b / a)**2
                t, i = _rootNd(_1_0 / r, 0, u, 0, w, g)  # eps
                a = _over(x, t * r + _1_0)
                b = _over(y, t     + _1_0)
                d =  hypot(x - a, y - b)
            else:  # on the ellipse
                a, b, d = x, y, _0_0
        else:  # x == 0
            if y < 0:
                b = -b
            a = x  # _copysign_0_0
            d = fabs(y - b)

    elif x:  # y == 0
        d, r = None, _over0(a * x, (a + b) * (a - b))
        if r:
            a *= r
            r = _1_0 - r**2
            if r > EPS02:
                b *= sqrt(r)
                d  = hypot(x - a, y - b)
        elif x < 0:
            a = -a
        if d is None:
            b = y  # _copysign_0_0
            d = fabs(x - a)

    else:  # x == y == 0
        a = x  # _copysign_0_0
        d = b

    return a, b, d, i


def _plumbTo5(x, y, z, Tun, eps=EPS):  # in .testTriaxials
    '''(INTERNAL) Nearest point on and distance to an I{un-/ordered} triaxial.

       @see: I{Eberly}'s U{Distance from a Point to ... an Ellipsoid ...<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    a, b, c, T = Tun._ordered4
    if Tun is not T:  # T is ordered, Tun isn't
        t = T._order3(x, y, z) + (T,)
        a, b, c, d, i = _plumbTo5(*t, eps=eps)
        return T._order3(a, b, c, reverse=True) + (d, i)

    if not (c > 0 and _isfinite(a)):
        raise _ValueError(a=a, b=b, c=c)

    if eps > 0:
        val = max(eps * 1e8, EPS)
    else:  # no validation
        val, eps = 0, max(EPS0, -eps)

    i = None
    if z:
        if y:
            if x:
                u =  fabs(x / a)
                v =  fabs(y / b)
                w =  fabs(z / c)
                g = _hypot2_1(u, v, w)
                if fabs(g) > EPS02:
                    r = T._c2_a2  # (c / a)**2
                    s = T._c2_b2  # (c / b)**2
                    t, i = _rootNd(_1_0 / r, _1_0 / s, u, v, w, g)  # eps
                    a = _over(x, t * r + _1_0)
                    b = _over(y, t * s + _1_0)
                    c = _over(z, t     + _1_0)
                    d =  hypot_(x - a, y - b, z - c)
                else:  # on the ellipsoid
                    a, b, c, d = x, y, z, _0_0
            else:  # x == 0
                a          =  x  # = _copysign_0_0(x)
                b, c, d, i = _plumbTo4(y, z, b, c, eps=eps)
        elif x:  # y == 0
            b          =  y  # = _copysign_0_0(y)
            a, c, d, i = _plumbTo4(x, z, a, c, eps=eps)
        else:  # x == y == 0
            if z < 0:
                c = -c
            a, b, d = x, y, fabs(z - c)

    else:  # z == 0
        u = _over0(a * x, T._a2c2)  # (a + c) * (a - c)
        v = _over0(b * y, T._b2c2)  # (b + c) * (b - c)
        s = _hypot2_1(u, v)
        if u and v and s < 0:
            a *= u
            b *= v
            c *= sqrt(-s)
            d  = hypot_(x - a, y - b, c)
        else:
            c          =  z  # _copysign_0_0(z)
            a, b, d, i = _plumbTo4(x, y, a, b, eps=eps)

    if val > 0:
        _validate(a, b, c, d, T, x, y, z, val)
    return a, b, c, d, i


def _reverseLatLon3(s, atan2_, v, forward_):
    '''(INTERNAL) Helper for C{.reverseBetaOmega} and C{.reverseLatLon}.
    '''
    x, y, z = s.xyz3
    d = hypot( x, y)
    a = atan2_(z, d)
    b = atan2_(y, x)
    h = v.minus_(*forward_(z, d, y, x)).length
    return a, b, (h or INT0)


def _rootNd(r, s, u, v, w, g, eps=EPS0):
    '''(INTERNAL) Robust 2-D or 3-D root finder: 2-D if C{s == v == 0} else 3-D root.

       @see: I{Eberly}'s U{Robust Root Finders ... and Listing 4<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    u *=  r
    v *=  s  # 0 for 2-D root
    t0 =  w - _1_0
    t1 = _0_0 if g < 0 else (hypot_(u, w, v) - _1_0)
    # assert t0 <= t1
    for i in range(1, _TRIPS):  # 48..58
        t = (t1 + t0) * _0_5
        e =  t1 - t0
        if eps > e > -eps or _isin(t, t0, t1):
            break
        g = fsumf_(_N_1_0,  # ~= _hypot2_1
                   _over02(u, t +  r),
                   _over02(w, t + _1_0), (
                   _over02(v, t +  s) if v else _0_0))
        if g > 0:
            t0 = t
        elif g < 0:
            t1 = t
        else:
            break
    else:  # PYCHOK no cover
        t = Fmt.no_convergence(e, eps)
        raise _ValueError(t, txt__=_rootNd)
    return t, i


def _tri_biaxial(tri_biax, where):
    '''(INTERNAL) Get a triaxail for C{tri_biax}.
    '''
    if isinstance(tri_biax, _UnOrderedTriaxialBase):
        T = tri_biax
    else:
        D = tri_biax if isinstance(tri_biax, Datum) else \
                  _spherical_datum(tri_biax, name__=where)  # typename
        T = D.ellipsoid._triaxial
    return T


def _validate(a, b, c, d, T, x, y, z, val):
    '''(INTERNAL) Validate an C{_plumTo5} result.
    '''
    e = T.sideOf(a, b, c, eps=val)
    if e:  # not near the ellipsoid's surface
        raise _ValueError(a=a, b=b, c=c, d=d,
                          sideOf=e, eps=val)
    if d:  # angle of delta and normal vector
        m = Vector3d(x, y, z).minus_(a, b, c)
        if m.euclid > val:
            m = m.unit()
            n = T.normal3d(a, b, c)
            e = n.dot(m)  # n.negate().dot(m)
            if not isnear1(fabs(e), eps1=val):
                raise _ValueError(n=n, m=m,
                                dot=e, eps=val)


if __name__ == _DMAIN_:

    from pygeodesy import printf
    from pygeodesy.interns import _COMMA_, _NL_, _NLATvar_

    T = Triaxial_(6378388.0, 6378318.0, 6356911.9461)
    t = T.height4(3909863.9271, 3909778.123, 3170932.5016)
    printf('# Bektas: %r', t)

    # __doc__ of this file, force all into registery
    t = [NN] + Triaxials.toRepr(all=True, asorted=True).split(_NL_)
    printf(_NLATvar_.join(i.strip(_COMMA_) for i in t))

# % python3 -m pygeodesy.triaxials
#
# Bektas: height4(x=3909251.554667, y=3909165.750567, z=3170432.501602, h=999.999996)

# **) MIT License
#
# Copyright (C) 2022-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
