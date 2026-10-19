
# -*- coding: utf-8 -*-

u'''I{Ordered} triaxial ellipsoid classes L{Triaxial3} and L{Triaxial3B} for conversion between
variuos lat-/longitudal and cartesian coordinates on a triaxial ellipsoid using L{Ang}, L{Deg},
L{Rad} lat-, longitude and heading angles.

Transcoded to pure Python from I{Karney}'s GeographicLib 2.7 C++ classes U{Ellipsoidal3<https://
GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Ellipsoidal3.html>} and U{Cartesian3
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Cartesian3.html>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2024-2025) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib 2.7 <https://GeographicLib.SourceForge.io/>}
documentation.

@var Triaxial3s.Amalthea: Triaxial3(name='Amalthea', a=125000, b=73000, c=64000, k2=0.106947697, kp2=0.893052303, volume=2446253479595252, area=93239507787.490356445, R2=86138.05359954)
@var Triaxial3s.Ariel: Triaxial3(name='Ariel', a=581100, b=577900, c=577700, k2=0.05866109, kp2=0.94133891, volume=812633172614203904, area=4211301462766.580078125, R2=578899.578791275)
@var Triaxial3s.Earth: Triaxial3(name='Earth', a=6378173.435, b=6378103.9, c=6356754.399999999, k2=0.996748146, kp2=0.003251854, volume=1083208241574987694080, area=510065911057440.9375, R2=6371008.987886564)
@var Triaxial3s.Enceladus: Triaxial3(name='Enceladus', a=256600, b=251400, c=248300, k2=0.369647336, kp2=0.630352664, volume=67094551514082248, area=798618496278.596679688, R2=252095.300756832)
@var Triaxial3s.Europa: Triaxial3(name='Europa', a=1564130, b=1561230, c=1560930, k2=0.093663002, kp2=0.906336998, volume=15966575194402123776, area=30663773697323.51953125, R2=1562096.533153486)
@var Triaxial3s.Io: Triaxial3(name='Io', a=1829400, b=1819300, c=1815700, k2=0.262045618, kp2=0.737954382, volume=25313121117889765376, area=41691875849096.734375, R2=1821464.812747882)
@var Triaxial3s.Mars: Triaxial3(name='Mars', a=3394600, b=3393300, c=3376300, k2=0.92878339, kp2=0.07121661, volume=162907283585817247744, area=144249140795107.4375, R2=3388064.624110653)
@var Triaxial3s.Mimas: Triaxial3(name='Mimas', a=207400, b=196800, c=190600, k2=0.359218713, kp2=0.640781287, volume=32587072869017956, area=493855762247.691833496, R2=198241.75359411)
@var Triaxial3s.Miranda: Triaxial3(name='Miranda', a=240400, b=234200, c=232900, k2=0.171062751, kp2=0.828937249, volume=54926187094835456, area=698880863325.757080078, R2=235828.692095158)
@var Triaxial3s.Moon: Triaxial3(name='Moon', a=1735550, b=1735324, c=1734898, k2=0.653331685, kp2=0.346668315, volume=21886698675223740416, area=37838824729886.09375, R2=1735257.329122863)
@var Triaxial3s.Tethys: Triaxial3(name='Tethys', a=535600, b=528200, c=525800, k2=0.243190549, kp2=0.756809451, volume=623086233855821440, area=3528073490771.393554688, R2=529863.348254881)
@var Triaxial3s.WGS84_3: Triaxial3(name='WGS84_3', a=6378171.36, b=6378101.609999999, c=6356751.84, k2=0.996738165, kp2=0.003261835, volume=1083207064030173855744, area=510065541435967.5, R2=6371006.679496506)
@var Triaxial3s.WGS84_35: Triaxial3(name='WGS84_35', a=6378172, b=6378102, c=6356752.314245179, k2=0.996726499, kp2=0.003273501, volume=1083207319768789942272, area=510065621722018.25, R2=6371007.180905545)
@var Triaxial3s.WGS84_3r: Triaxial3(name='WGS84_3r', a=6378172, b=6378102, c=6356752, k2=0.996726547, kp2=0.003273453, volume=1083207266220584468480, area=510065604942135.875, R2=6371007.076110449)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.angles import Ang, Ang_, _Ang3Tuple,  atan2, sincos2, _SinCos2
from pygeodesy.basics import _copysign, map1
from pygeodesy.constants import EPS, EPS_2, EPS02, EPS8, INT0, NAN, \
                               _EPSqrt, _SQRT3, _copysign_0_0, _copysign_1_0, \
                               _flipsign, _isfinite, _over, _1_over, _0_0, \
                               _0_5, _N_1_0, _1_0, _2_0, _3_0, _4_0, _9_0
from pygeodesy.errors import _xattr, _xkwds, _xkwds_get, _xkwds_pop2
from pygeodesy.fmath import cbrt2, fdot, hypot, hypot2, norm2,  fabs, sqrt
from pygeodesy.fsums import Fsum, fsumf_,  Fmt
from pygeodesy.interns import NN, _azimuth_, _DMAIN_, _h_, _lam_, _lat_, \
                             _lon_, _phi_
# from pygeodesy.lazily import _ALL_LAZY  # from .vector3d
# from pygeodesy.named import _Pass  # from .namedTuples
from pygeodesy.namedTuples import _NamedTuple, _xyzh_Tuple,  _Pass, Property_RO
# from pygeodesy.props import Property_RO  # from .namedTuples
# from pygeodesy.streprs import Fmt  # from .fsums
from pygeodesy.triaxials.bases import _bet_, _HeightINT0, LLK, _llk_, \
                                      _MAXIT, _omg_, _otherV3d_, _sqrt0, \
                                      _Triaxial3Base, TriaxialError, \
                                      _TriaxialsBase
from pygeodesy.units import Degrees, Lat, Lon, Meter, Radians, Radius_
# from pygeodesy.utily import atan2, sincos2  # from .triaxials.angles
from pygeodesy.vector3d import Vector3d,  _ALL_LAZY

# from math import fabs, sqrt  # from .fmath
from random import random

__all__ = _ALL_LAZY.triaxials_triaxial3
__version__ = '26.08.06'

_alp_  = 'alp'
_NAN3d =  Vector3d(NAN, NAN, NAN)
_TOL   =  cbrt2(EPS)
_TOL2  = _TOL**2  # cbrt(EPS)**4
_zet_  = 'zet'
_27_0  =  27.0


class BetOmgAlp5Tuple(_Ang3Tuple):
    '''5-Tuple C{(bet, omg, alp, h, llk)} with I{ellipsoidal}
       lat- C{bet}, longitude C{omg} and azimuth C{alp}, all
       in L{Ang}les on and height C{h} off the triaxial's
       surface and kind C{llk} set to C{LLK.ELLIPSOIDAL}.
    '''
    _Names_ = (_bet_, _omg_, _alp_, _h_,         _llk_)
    _Units_ = ( Ang,   Ang,  _Pass, _HeightINT0, _Pass)


class Cartesian5Tuple(_xyzh_Tuple):
    '''5-Tuple C{(x, y, z, h, llk)} with I{cartesian} C{x}, C{y}
       and C{z} coordinates on and height C{h} above or below the
       triaxial's surface and kind C{llk} set to the original
       C{LLK} or C{None}.
    '''
    _Names_ = _xyzh_Tuple._Names_ +            (_llk_,)
    _Units_ = (Meter, Meter, Meter, _HeightINT0, _Pass)

    def __new__(cls, x, y, z, h=0, llk=None, **kwds):  # **iteration_name
        args = x, y, z, (h or INT0), llk
        return _xyzh_Tuple.__new__(cls, args, **kwds)


class _Fp2(object):
    '''(INTERNAL) Function and derivate evaluation.
    '''
    def __init__(self, rs, ls, n=1):
        # assert 0 < n <= 2
        self._2   = n == 2
        self._rls = tuple((p, q) for p, q in zip(rs, ls) if p)

    def __call__(self, p):
        # Evaluate C{f(p) = sum((rs[k] / (p + ls[k]))**n,
        # k=0..2) - 1} and its derivative C{fp}.
        f  = _N_1_0
        fc =  fp = _0_0
        _D =  EPS_2
        _2 =  self._2
        for g, q in self._rls:
            q  = _1_over(p + q)
            g *=  q
            if _2:
                g *= g
                q += q
            r   = round(g / _D) * _D
            f  += r
            fc += g - r
            fp -= g * q
        return (f + fc), fp


class LatLonAzi5Tuple(_NamedTuple):
    '''5-Tuple C{(lat, lon, azimuth, h, llk)} with triaxial
       C{lat}-, C{lon}gitude and C{azimuth}, all in C{degrees},
       height C{h} off the triaxial's surface and kind C{llk}
       either C{LLK.GEODETIC} or C{LLK.GEODESIC_LON0}.
    '''
    _Names_ = (_lat_, _lon_, _azimuth_, _h_,         _llk_)
    _Units_ = ( Lat,   Lon,  _Pass,     _HeightINT0, _Pass)


class PhiLamZet5Tuple(_Ang3Tuple):
    '''5-Tuple C{(phi, lam, zet, h, llk)} with triaxial lat-
       lat- C{phi}, longitude C{lam} and azimuth C{zet}, all
       in L{Ang}les on and height C{h} off the triaxial's
       surface and kind C{llk} set to an C{LLK}.
    '''
    _Names_ = (_phi_, _lam_, _zet_, _h_,         _llk_)
    _Units_ = ( Ang,   Ang,  _Pass, _HeightINT0, _Pass)

    def toLatLonAzi5(self, **name):
        '''Return this tuple as L{LatLonAzi5Tuple} with
           C{lat}, C{lon} and C{azimuth} in C{degrees} or
           C{azimuth} C{INT0}.
        '''
        u =  self  # .toUnits(Error=TriaxialError)
        n = _xkwds_get(name, name=u.name)
        z =  u.zet  # None or isAng with turns
        z =  INT0 if z is None else z.degrees
        t =  u.phi.degrees0, u.lam.degrees0, z, u.h, u.llk
        return LatLonAzi5Tuple(t, name=n or NN)


class Triaxial3(_Triaxial3Base):
    '''I{Ordered} triaxial ellipsoid convering between cartesian and
       lat-/longitudes using using class L{Ang}.

       @see: L{Triaxial<triaxials.triaxial5.Triaxial>} for details.
    '''
    def _cardinal2(self, v, mer, llk):  # cardinaldir
        '''(INTERNAL) Get 2-tuple C{(n, e)} at C{mer}idian.
        '''
        # assert isinstance(v, Vector3d) and isinstance(mer, Ang) \
        #                                and isinstance(llk, LLK.__class__)
        a2, b2, c2 = self._a2b2c23
        if llk._X:
            a2, c2 = c2, a2
            v = v._roty(True)  # +1
        x, y, z = v.xyz3
        if x or y:
            s = (-z) / c2
            z = x**2 / a2 + y**2 / b2
        else:
            y, x, _ = mer.scn3
            s = _copysign_1_0(-z)
            z = _0_0
        n = Vector3d(x * s, y * s, z).unit()
        e = v.dividedBy_(a2, b2, c2).unit()  # normvec
        e = n.cross(e).unit()
        if llk._X:
            e = e._roty(False)  # -1
            n = n._roty(False)  # -1
        return n, e

    def forward(self, lat, lon, **height_unit_name):
        '''Convert a I{geodetic} lat- and longitude to a cartesian
           on this triaxial's surface.

           @arg lat: Geodetic latitude (C{degrees} or B{C{unit}}).
           @arg lon: Geodetic longitude (C{degrees} or B{C{unit}}).
           @kwarg height_unit_name: Optional C{B{height}=0} (C{meter}),
                         scalar C{B{unit}=}L{Radians} (or L{Degrees})
                         and C{B{name}=NN} (C{str}).

           @return: A L{Cartesian5Tuple}C{(x, y, z, h, llk)} with
                    C{h=B{height}} and kind C{llk=LLK.GEODETIC} or
                    C{llk=LLK.GEODETIC_LON0}.

           @see: Methods L{Triaxial3.forwardPhiLam} and L{Triaxial3.reverse}.
        '''
        llk = LLK.GEODETIC_LON0 if self.Lon0 else LLK.GEODETIC
        height_unit_name = _xkwds(height_unit_name, unit=Degrees)
        return self.forwardPhiLam(lat, lon, llk=llk, **height_unit_name)

    def forwardBetOmg(self, bet, omg, height=0, **unit_name):  # elliptocart2
        '''Convert an I{ellipsoidal} lat- and longitude to a cartesian
           on this triaxial's surface.

           @arg bet: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omg: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg height: Height above or below this triaxial's surface (C{meter},
                          same units as this triaxial's semi-axes).
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar
                             C{B{unit}=}L{Radians} (or L{Degrees}).

           @return: A L{Cartesian5Tuple}C{(x, y, z, h, llk)} with C{h=B{height}}
                    and kind C{llk=LLK.ELLIPSOIDAL}.

           @see: Method L{Triaxial3.reverseBetOmg}.
        '''
        ct, _ = self.forwardBetOmgAlp2(bet, omg, None, height, **unit_name)
        return ct

    forwardBetaOmega = forwardBetOmg  # for backward compatibility

    def forwardBetaOmega_(self, sbeta, cbeta, somega, comega, **name):
        '''DEPRECATED on 2025.11.15, like C{Triaxial.forwardBetaOmega_}.'''
        return self.forwardBetaOmega(Ang_(sbeta, cbeta),
                                     Ang_(somega, comega), **name)

    def forwardBetOmgAlp2(self, bet, omg, alp, height=0, **unit_name):  # elliptocart2
        '''Convert an I{ellipsoidal} lat-, longitude and heading to a
           cartesian and a direction on this triaxial's surface.

           @arg bet: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omg: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @arg alp: Azimuth of the heading (C{Ang}, B{C{unit}} or C{None}).
           @kwarg height: Height above or below this triaxial's surface (C{meter},
                          same units as this triaxial's semi-axes).
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}), scalar
                       C{B{unit}=}L{Radians} (or L{Degrees}).

           @return: 2-Tuple C{(cartesian, direction)} with C{cartesian} a
                    L{Cartesian5Tuple}C{(x, y, z, h, llk)} with C{h=B{height}},
                    kind C{llk=LLK.ELLIPSOIDAL} and C{direction} a C{Vector3d}
                    tangent to this triaxial's surface or C{None}.

           @see: Method L{Triaxial3.reverseBetOmgAlp}.
        '''
        h, llk, unit, name = _h_llk_unit_name(height, **unit_name)
        a, b, c = self._abc3
        if h:  # Cartesian.elliptocart
            a, b, _ = self._a2b2c23
            h = _HeightINT0(h)
            s = (c * _2_0 + h) * h
            if s < 0:
                s = -min(a, b, -s)
            a  = sqrt(a + s)
            b  = sqrt(b + s)
            c += h
        sb, cb = _SinCos2(bet, unit)
        so, co = _SinCos2(omg, unit)
        k,  kp =  self._k_kp
        tx, tz = _txtz2(cb, so, k, kp)
        ct = Cartesian5Tuple(a * co * tx,
                             b * cb * so,
                             c * sb * tz,
                             h, llk, **name)

        if alp is None:  # or h?
            dir3d = None  # _NAN3d?
        else:
            try:
                sa, ca = _SinCos2(alp, unit)
            except Exception as X:
                raise TriaxialError(alp=alp, cause=X)
            a, b, c =  self._abc3
            if k and kp and not (cb or so):
                c2s2_b = (ca - sa) * (ca + sa) / b
                dir3d = Vector3d(a  *  k * co * c2s2_b
                                -co * sb * ca * sa * _2_0,
                                 c  * kp * sb * c2s2_b)
            else:
                if not tx:  # at oblate pole tx -> |cos(bet)|
                    c = _flipsign(co, cb)
                    n =  Vector3d(-c  * sb,
                                  -so * sb, _0_0)
                    e =  Vector3d(-so,   c, _0_0)
                elif not tz:  # at prolate pole tz -> |sin(omg)|
                    s = _flipsign(sb, so)
                    n =  Vector3d(_0_0,      -s, cb)
                    e =  Vector3d(_0_0, cb * co, co * s)
                else:
                    k2, kp2 = self._k2_kp2
                    n = Vector3d(-a * k2  * sb * cb * co / tx,
                                 -b * sb  * so,   c * cb * tz)
                    e = Vector3d(-a * tx  * so,   b * cb * co,
                                  c * kp2 * sb * so * co / tz)
                dir3d  = n.unit().times(ca)  # NAN
                dir3d += e.unit().times(sa)  # NAN
            dir3d.name = ct.name
        return ct, dir3d

    def forwardCartesian(self, x_ct, y=None, z=None, normal=True, **eps_llk_name):
        '''Project any cartesian I{onto} this triaxial's surface.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg normal: If C{True}, the projection is C{perpendicular} to the surface,
                          otherwise C{radial} to the center of this triaxial (C{bool}).
           @kwarg eps_llk_name: Root finder tolerance C{B{eps}=EPS}, kind C{B{llk}=None}
                      overriding C{B{x_ct}.llk} and optional C{B{name}="height4"} (C{str}).

           @return: A L{Cartesian5Tuple}C{(x, y, z, h, llk)}.

           @see: Method L{Triaxial3.reverseCartesian} to reverse the projection and
                 function L{height4<triaxials.triaxial5.height4>} for more details.
        '''
        llk, kwds = _xkwds_pop2(eps_llk_name, llk=_xattr(x_ct, llk=None))
        h = self.sideOf(x_ct, y, z)
        if h:  # signed, square
            v = self.height4(x_ct, y, z, normal=normal, **kwds)
            h = v.h
        else:  # on the surface
            v = _otherV3d_(x_ct, y, z)
        n = _xkwds_get(kwds, name=NN)
        return Cartesian5Tuple(v.x, v.y, v.z, h, llk, iteration=v.iteration, name=n)

    def forwardLatLon(self, lat, lon, height=0, llk=LLK.ELLIPSOIDAL, **unit_name):  # anytocart2
        '''Convert any lat-/longitude kind to a cartesian on this triaxial's surface.

           @arg lat: Latitude (C{Ang} or B{C{unit}}).
           @arg lon: Longitude (C{Ang} or B{C{unit}}).
           @kwarg height: Height above or below this triaxial's surface (C{meter}, same
                          units as this triaxial's semi-axes).
           @kwarg llk: The kind (an L{LLK}).
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar C{B{unit}=}L{Degrees}
                             (or L{Radians}).

           @return: A L{Cartesian5Tuple}C{(x, y, z, h, llk)} with height C{h=B{height}} and
                    kind C{llk=B{llk}}.

           @see: Method L{Triaxial3.reverseLatLon}.
        '''
        _fwd = self.forwardBetOmg if llk in LLK._NOIDAL else \
               self.forwardPhiLam  # PYCHOK OK
        return _fwd(lat, lon, height=height, llk=llk, **_xkwds(unit_name, unit=Degrees))

    def forwardPhiLam(self, phi, lam, height=0, llk=LLK.GEODETIC, **unit_name):  # generictocart2
        '''Convert any lat-/longitude kind to a cartesian on this triaxial's surface.

           @arg phi: Latitude (C{Ang} or B{C{unit}}).
           @arg lam: Longitude (C{Ang} or B{C{unit}}).
           @kwarg height: Height above or below this triaxial's surface (C{meter}, same
                          units as this triaxial's semi-axes).
           @kwarg llk: The kind (an L{LLK}).
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar C{B{unit}=}L{Radians}
                             (or L{Degrees}).

           @return: A L{Cartesian5Tuple}C{(x, y, z, h, llk)} with height C{h=B{height}} and
                    kind C{llk=B{llk}}.

           @note: Longitude C{B{lam} -= Lon0} if C{B{llk} is LLK.GEODETIC_LON0}.

           @see: Method L{Triaxial3.reverseLatLon}.
        '''
        ct, _ = self.forwardPhiLamZet2(phi, lam, None, height=height, llk=llk, **unit_name)
        return ct

    def forwardPhiLamZet2(self, phi, lam, zet, height=0, llk=LLK.GEODETIC, **unit_name):  # generictocart2
        '''Convert a lat-, longitude and heading to a cartesian and a direction
           on this  triaxial's surface.

           @arg phi: Latitude (C{Ang} or B{C{unit}}).
           @arg lam: Longitude (C{Ang} or B{C{unit}}).
           @arg zet: Azimuth of the heading (C{Ang}, B{C{unit}} or C{None}).
           @kwarg height: Height above or below this triaxial's surface (C{meter},
                          same units as this triaxial's semi-axes).
           @kwarg llk: The kind (an L{LLK}).
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar
                       C{B{unit}=}L{Radians} (or L{Degrees}).

           @return: 2-Tuple C{(cartesian, direction)} with the C{cartesian} a
                    L{Cartesian5Tuple}C{(x, y, z, h, llk)} with height C{h=0},
                    kind C{llk=B{llk}} and C{direction}, a C{Vector3d} on and
                    tangent to this triaxial's surface.

           @note: Longitude C{B{lam} -= Lon0} if C{B{llk} is LLK.GEODETIC_LON0}.

           @see: Method L{Triaxial3.reversePhiLamZet}.
        '''
        unit, name = _xkwds_pop2(unit_name, unit=Radians)
        try:
            sa, ca = _SinCos2(phi, unit)
            if llk is LLK.GEODETIC_LON0 and self.Lon0:
                lam  = Ang.fromScalar(lam, unit=unit)
                lam -= self.Lon0
            sb, cb = _SinCos2(lam, unit)
        except Exception as X:
            raise TriaxialError(phi=phi, lam=lam, llk=llk, cause=X)
        v, _, llk, name = _v_h_llk_name(ca * cb, ca * sb, sa, llk=llk, **name)
        if llk and llk._X:
            v = v._roty(False)  # -1
        d, t = _d_t(self, llk)
        if t:
            v = v.times_(*t)
        if d:
            d = v.dividedBy_(*self._abc3).length
            v = v.dividedBy(d)

        h = _HeightINT0(height)
        if h:  # cart2cart
            v, h = self._toHeight2(v, h)
        ct = Cartesian5Tuple(v.x, v.y, v.z, h, llk, **name)

        if zet is None:
            dir3d = None
        else:
            try:
                s, c = _SinCos2(zet, unit)
            except Exception as X:
                raise TriaxialError(zet=zet, cause=X)
            n, e   = self._meridian2(v, lam, llk)
            dir3d  = n.times(c)
            dir3d += e.times(s)
            dir3d.name = ct.name
        return ct, dir3d

    def _meridian(self, lam, llk):
        '''(INTERNAL) Get the meridian plane's at C{lam}.
        '''
        _, t = _d_t(self, llk)
        if t:
            a, b, c = t
            lam = lam.mod((c if llk._X else a) / b)
        return lam

    def _meridian2(self, v, lam, llk):
        '''(INTERNAL) Get 2-tuple C{(n, e)} at C{lam} meridian.
        '''
        mer =  self._meridian(lam, llk)
        return self._cardinal2(v, mer, llk)

    def normed2(self, x_ct, y=None, z=None, dir3d=None, **llk_name):  # Ellipsoid3.Norm
        '''Scale a cartesian and direction to this triaxial's surface.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg dir3d: The direction (C{Vector3d} or C{None}).
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                      overriding C{B{x_ct}.llk}.

           @return: 2-Tuple C{(cartesian, direction)} with the C{cartesian} a
                    L{Cartesian5Tuple}C{(x, y, z, h, llk)} and C{direction}, a
                    C{Vector3d} tangent to this triaxial's surface or C{None}
                    iff C{B{dir3d} is None}.
        '''
        v, h, llk, name = _v_h_llk_name(x_ct, y, z, **llk_name)

        u  = v.dividedBy_(*self._abc3).length
        r  = v.dividedBy(u) if u else _NAN3d
        ct = Cartesian5Tuple(r.x, r.y, r.z, h, llk, **name)

        if isinstance(dir3d, Vector3d):
            if u:  # and r is not _NAN3d
                u = r.dividedBy_(*self._a2b2c23)
                d = dir3d.dot(u)
                if _isfinite(d) and u.length2:
                    u = u.times(d / u.length2)
                    dir3d = dir3d.minus(u).unit()  # NAN
                else:
                    dir3d = _NAN3d
            else:
                dir3d = _NAN3d
            dir3d.name = ct.name
        return ct, dir3d

    def reverse(self, x_ct, y=None, z=None, **name):
        '''Convert a cartesian I{on this triaxial's surface} to I{geodetic} lat-
           and longitude.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{LatLonAzi5Tuple}C{(lat, lon, azimuth, h, llk)} with C{azimuth}
                    usually C{INT0}, height C{h} off this triaxial's surface in C{meter}
                    or C{INT0} and C{llk} either C{LLK.GEODETIC} or C{LLK.GEODETIC_LON0}.
        '''
        llk = LLK.GEODETIC_LON0 if self.Lon0 else LLK.GEODETIC
        r = self.reversePhiLamZet(x_ct, y, z, llk=llk, **name)
        return r.toLatLonAzi5(name=r.name)

    def reverseBetOmg(self, x_ct, y=None, z=None, **llk_name):  # Cartesian3.carttoellip
        '''Convert a cartesian I{on this triaxial's surface} to an I{ellipsoidal}
           lat-/longitude.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                      overriding C{B{x_ct}.llk}.

           @return: A L{BetOmgAlp5Tuple}C{(bet, omg, alp, h, llk)} with C{alp=None}
                    and C{llk=LLK.ELLIPSOIDAL}.
        '''
        v, _, llk, name = _v_h_llk_name_NOIDAL(x_ct, y, z, **llk_name)

        _, y2, z2 = rs = v.x2y2z23
        l0, l1, _ = ls = self._lcc23
        qmax = fsumf_(*rs)
        qmin = q = max(z2, y2 + z2 - l1, qmax - l0)
        _fp2 = _Fp2(rs, ls, n=1)
        f, _ = _fp2(q)
        if f > _TOL2:  # neg means convergence
            q = max(qmin, min(qmax, _cubic(rs, qmax, l0, l1)))
            f, fp = _fp2(q)
            if fabs(f) > _TOL2:
                q =  max(qmin, q - _over(f, fp))
                q = _solve(_fp2, q, self.b2)

        a, b, c = map1(_sqrt0, l0 + q, l1 + q, q)  # axes (a, b, c)
        h = (c - self.c) or INT0
        bet, omg, _ = self._reverseBetOmgAlp3(v, None, a, b, c, **name)
        return BetOmgAlp5Tuple(bet, omg, None, h, llk, **name)

    reverseBetaOmega = reverseBetOmg  # for backward compatibility

    def reverseBetOmgAlp(self, x_ct, y=None, z=None, dir3d=None, **llk_name):  # Ellipsoid3.cart2toellip[int]
        '''Convert a cartesian and direction I{on this triaxial's surface} to an
           I{ellipsoidal} lat-, longitude and heading.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_ct}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg dir3d: The direction (C{Vector3d} or C{None}).
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                      overriding C{B{x_ct}.llk}.

           @return: A L{BetOmgAlp5Tuple}C{(bet, omg, alp, h, llk)} with C{alp=None}
                    if C{B{dir3d} is None} and C{llk=LLK.ELLIPSOIDAL}.
        '''
        v, h, llk, name = _v_h_llk_name_NOIDAL(x_ct, y, z, **llk_name)
        bet, omg, alp = self._reverseBetOmgAlp3(v, dir3d, **name)
        return BetOmgAlp5Tuple(bet, omg, alp, h, llk, **name)

    def _reverseBetOmgAlp3(self, v, dir3d, *a_b_c, **name):  # cart2toellipint
        '''(INTERNAL) Helper for methods C{reverseBetOmg/-Alp}.
        '''
        k,  kp  = self._k_kp
        k2, kp2 = self._k2_kp2
        a, b, c = a_b_c or self._abc3
        V       = v.dividedBy_(a, b, c)
        X, E, Z = V.xyz3  # Xi, Eta, Zeta
        h = fabs(E * k * kp * _2_0)
        if v.y or fabs(v.x) != a * kp2 or \
                  fabs(v.z) != c * k2:
            g = fdot(V.x2y2z23, k2, (k2 - kp2), -kp2)
            h = hypot(g, h)
        else:
            g = _0_0
        if h < EPS02:
            so = cb = _0_0
        elif g < 0:
            h  = _over(sqrt((h - g) * _0_5), kp)
            so = _copysign(h, E)
            cb =  fabs(_over(E, so))
        else:
            cb = _over(sqrt((h + g) * _0_5), k)
            so = _over(E, cb)
        tx, tz = _txtz2(cb, so, k, kp)
        sb = (Z / tz) if tz else _N_1_0
        co = (X / tx) if tx else _1_0
        bet = Ang_(sb, cb, **name)
        omg = Ang_(so, co, **name)

        if isinstance(dir3d, Vector3d):  # cart2toellip(bet, omg, V) -> alp
            if cb or so or not (tx and tz):  # not umbilical
                if not tx:
                    n = Vector3d(-co, -so, tx) * sb
                    e = Vector3d(-so,  co, tx)
                elif not tz:
                    n = Vector3d(tz, -sb, cb)
                    e = Vector3d(tz,  cb, sb) * co
                else:
                    n = Vector3d(-a * sb * k2 * cb * co / tx,
                                 -b * sb * so,   c * cb * tz)
                    e = Vector3d(-a * so * tx,    b * cb * co,
                                  c * so * kp2 * sb * co / tz)
                sa = dir3d.dot(e.unit())  # NAN
                ca = dir3d.dot(n.unit())  # NAN
            else:  # at umbilicial PYCHOK no cover
                x, z = norm2(co * tx / a, sb * tz / c)  # _MODS.karney._norm2
                v   =  dir3d * (sb * co)  # dir3d.times(sb * co)
                s2a = -v.y
                c2a =  fdot(v, z, 0, -x)  # v.x * z - v.z * x
                sa  =  ca = -sb
                sa *= _copysign(_1_0 - c2a, s2a) if c2a < 0 else  s2a
                ca *=  fabs(s2a)                 if c2a < 0 else (c2a + _1_0)
            alp = Ang_(sa, ca, **name)
        elif dir3d is None:
            alp = None  # Ang.NAN(**name)
        else:
            raise TriaxialError(dir3d=dir3d)
        return bet, omg, alp

    def reverseCartesian(self, x_ct, y=None, z=None, height=0, normal=True, **llk_name):  # cart2tocart
        '''"Unproject" a cartesian I{off} this triaxial's surface.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg height: Height above or below this triaxial's surface (C{meter},
                          same units as this triaxial's semi-axes).
           @kwarg normal: If C{True}, B{C{height}} is C{perpendicular} to the surface,
                          otherwise C{radial} to the center of this triaxial (C{bool}).
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}}
                            overriding C{B{x_ct}.llk}.

           @return: L{Cartesian5Tuple}C{(x, y, z, h, llk)}.

           @raise TrialError: Cartesian B{C{x_ct}} or C{(x, y, z)} not on this
                              triaxial's surface.

           @see: Methods L{Triaxial3.forwardCartesian}.
        '''
        kwds = _xkwds(llk_name, llk=_xattr(x_ct, llk=None))
        v, _, llk, name = _v_h_llk_name(x_ct, y, z, **kwds)
        _ =  self._sideOn(v)
        h = _HeightINT0(height)
        if h:
            v, h = self._toHeight2(v, h, normal)
        return Cartesian5Tuple(v.x, v.y, v.z, h, llk, **name)

    def reverseLatLon(self, x_ct, y=None, z=None, **llk_name):  # cart2toany
        '''Convert a cartesian I{on this triaxial's surface} to a lat-/longitude.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                            overriding C{B{x_ct}.llk}.

           @return: A L{BetOmgAlp5Tuple}C{(bet, omg, alp, h, llk)} with C{alp=None} or
                    a L{PhiLamZet5Tuple}C{(phi, lam, zet, h, llk)} with C{zet=None}.

           @note: Longitude C{B{lam} += Lon0} if C{B{llk} is LLK.GEODETIC_LON0}.
        '''
        llk, name = _xkwds_pop2(llk_name, llk=_xattr(x_ct,
                                          llk=LLK.ELLIPSOIDAL))
        _rev = self.reverseBetOmg if llk in LLK._NOIDAL else \
               self.reversePhiLam  # PYCHOK OK
        return _rev(x_ct, y, z, llk=llk, **name)

    def reversePhiLam(self, x_ct, y=None, z=None, **llk_name):  # cart2togeneric
        '''Convert a cartesian I{on this triaxial's surface} to lat-/longitude.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                            overriding C{B{x_ct}.llk}.

           @return: A L{PhiLamZet5Tuple}C{(phi, lam, zet, h, llk)} with C{zet=None}.

           @note: Longitude C{B{lam} += Lon0} if C{B{llk} is LLK.GEODETIC_LON0}.
        '''
        return self.reversePhiLamZet(x_ct, y, z, **llk_name)

    def reversePhiLamZet(self, x_ct, y=None, z=None, dir3d=None, **llk_name):  # cart2togeneric(R, V, ...
        '''Convert a cartesian and direction to lat-, longitude and azimuth.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg dir3d: Optional direction (C{Vector3d} or C{None}).
           @kwarg llk_name: Optional C{B{name}=NN} (C{str}) and kind C{B{llk}=None}
                            overriding C{B{x_ct}.llk}.

           @return: A L{PhiLamZet5Tuple}C{(phi, lam, zet, h, llk)} with C{zet=None}
                    if C{B{dir3d} is None}.

           @note: Longitude C{B{lam} += Lon0} if C{B{llk} is LLK.GEODETIC_LON0}.
        '''
        ct = self.toTriaxial5(x_ct, y, z, h=NAN, **llk_name)
        v, h, llk, name = _v_h_llk_name(ct)
        _, t = _d_t(self, llk)
        if t:
            v = v.dividedBy_(*t)
        if llk._X:
            v = v._roty(True)  # +1
        phi = Ang_(v.z, hypot(v.x, v.y), **name)
        lam = Ang_(v.y, v.x, **name)  # Ang(0, 0) -> 0

        if dir3d is None:
            zet = None
        elif isinstance(dir3d, Vector3d):
            n, e = self._meridian2(v, lam, llk)
            zet  = Ang_(dir3d.dot(e),
                        dir3d.dot(n), **name)
        else:
            raise TriaxialError(dir3d=dir3d)
        if llk is LLK.GEODETIC_LON0 and self.Lon0:
            lam += self.Lon0
        return PhiLamZet5Tuple(phi, lam, zet, h, llk, **name)

    def random2(self, llk=LLK.ELLIPSOIDAL, both=False, _rand=random):
        '''Return a random cartesian with/out direction on this triaxial's surface.

           @kwarg llk: The kind (an L{LLK}).
           @kwarg both: If C{True}, generate a random direction (C{bool}).

           @return: 2-Tuple C{(cartesian, direction)} with the C{cartesian} a
                    L{Cartesian5Tuple}C{(x, y, z, h, llk)} and C{direction}, a
                    C{Vector3d} tangent to this triaxial's surface or C{None}
                    iff C{B{both} is False}.
        '''
        for _ in range(_MAXIT):
            for _ in range(_MAXIT):
                v = Vector3d(_rand(), _rand(), _rand())
                u = v.length
                if u and _isfinite(u):
                    break
            else:
                raise TriaxialError(Fmt.no_convergence(u))  # PYLINT E1120
            v = v.dividedBy(u).times_(*self._abc3)  # PYLINT E1101
            q = v.dividedBy_(*self._a2b2c23).length * self.c
            if 0 < q <= _1_0:  # _uni(q) < q:
                break
        else:
            raise TriaxialError(Fmt.no_convergence(q))  # PYLINT E1120
        ct = Cartesian5Tuple(v.x, v.y, v.z, INT0, llk, name__=self.random2)
        v  = None
        if both:
            for _ in range(_MAXIT):
                v = Vector3d(_rand(), _rand(), _rand())
                u = v.length
                if u:
                    u = v.dividedBy(u).dividedBy_(*self._a2b2c23)
                    d = v.dot(u) / u.length2
                    v = v.minus(u.times(d))
                    u = v.length  # normvec
                    if u and _isfinite(u):
                        v = v.dividedBy(u)
                        break
            else:
                raise TriaxialError(Fmt.no_convergence(u))  # PYLINT E1120
            v.name = ct.name
        return ct, v

    def _toHeight2(self, v, h, normal=True):
        '''(INTERNAL) Move cartesian C{Vector3d B{v}} to height C{h}.
        '''
        n = v.dividedBy_(*self._a2b2c23) if normal else v
        if n.length > EPS02:
            h = max(h, -self.c)
            v = v.plus(n.times(h / n.length))
        return v, h

    def toOther(self, lat, lon, llk1=LLK.GEODETIC, llk2=LLK.GEODETIC, **unit_name):  # anytoany
        '''Convert one lat-/longitude kind to an other.

           @arg lat: Latitude (C{Ang} or B{C{unit}}).
           @arg lon: Longitude (C{Ang} or B{C{unit}}).
           @kwarg llk1: The given kind (an L{LLK}).
           @kwarg llk2: The result kind (an L{LLK}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{BetOmgAlp5Tuple}C{(bet, omg, alp, h, llk)} with C{alp=None} or
                    a L{PhiLamZet5Tuple}C{(phi, lam, zet, h, llk)} with C{zet=None}.

           @see: Methods L{Triaxial3.forwardLatLon} and -L{reverseLatLon}.
        '''
        ct = self.forwardLatLon(lat, lon, llk=llk1, **unit_name)
        r  = self.reverseLatLon(ct,       llk=llk2, name=ct.name)
#       a, b = r[:2]
#       if not isAng(lat):
#           a = float(a)
#       if not isAng(lon):
#           b = float(b)
#       if (a, b) =! r[:2]:
#           r = r._dup(a, b)
        return r

    def toTriaxial5(self, x_ct, y=None, z=None, **triaxial_h_llk_name):  # carttocart2
        '''Find the closest cartesian on this or on another triaxial's surface.

           @arg x_ct: X component (C{scalar}) or a cartesian (L{Cartesian5Tuple} or
                      any C{Cartesian}, L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple}
                      or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg triaxial_llk_name: Optional C{B{triaxial}=self} (C{Triaxial3}),
                           C{B{name}=NN} (C{str}), height C{B{h}} and kind C{B{llk}}
                           overriding C{B{x_ct}.h} respectively C{B{x_ct}.llk}.

           @return: L{Cartesian5Tuple}C{(x, y, z, h, llk)}

           @raise TriaxialError: If C{B{triaxial}} is not a L{Triaxial3}.

           @see: Functions L{hartzell4<triaxials.triaxial5.hartzell4>} and
                 L{height4<triaxials.triaxial5.height4>} and methods.
        '''
        T, name = _xkwds_pop2(triaxial_h_llk_name, triaxial=self)
        if not isinstance(T, Triaxial3):
            raise TriaxialError(triaxial=T, x=x_ct, y=y, z=z)

        v, h, llk, name = _v_h_llk_name(x_ct, y, z, **name)
        if h or T is not self:
            l0, l1, _ = ls = T._lcc23
            r =  Vector3d(*T._ztol(v))
            s =  r.times_(*T._abc3)
            p =  max(fabs(s.z), hypot(s.x, s.y) - l1, s.length - l0)  # PYLINT E1101
            h = _solve(_Fp2(s.xyz3, ls, n=2), p, T.b2)  # PYLINT E1101
            v =  r.times_(*(_over(n, h + l_) for n, l_ in zip(T._a2b2c23, ls)))
            if not h:  # handle h == 0, v.y indeterminate
                x = v.x if l0 else r.x  # sphere
                y = v.y if l1 else r.y  # sphere or prolate
                s = _1_0 - hypot2(x / T.a, y / T.b)
                z = (sqrt(s) * _copysign(T.c, r.z)) if s > EPS02 else _0_0
                v = Vector3d(x, y, z)
            h -= T.c2
            if h and v.length:
                h *= v.dividedBy_(*T._a2b2c23).length
        return Cartesian5Tuple(v.x, v.y, v.z, (h or INT0), llk, **name)

    @Property_RO
    def _ZTOL(self):
        return self.b * (EPS / 8)

    def _ztol(self, v):
        for x in v.xyz3:
            yield x if fabs(x) > self._ZTOL else _copysign_0_0(x)


class Triaxial3B(Triaxial3):
    '''Triaxial ellipsoid specified by its middle semi-axis and shape.

       @see: L{Triaxial3} for more information.
    '''
    def __init__(self, b, e2=_0_0, k2=_1_0, kp2=_0_0, **name):  # PYLINT W0231
        '''New, L{Triaxial3B} instance.

           @see: L{Triaxial<triaxials.triaxial5.Triaxial>} for details.
        '''
        self._init_abc3_e2_k2_kp2(Radius_(b=b), e2, k2, kp2, **name)


def _cubic(rs, rt, l0, l1):  # Cartesian3.cubic
    '''(INTERNaL) Solve sum(R2[i]/(z + lq2[i]), i=0,1,2) - 1 = 0
        with lq2[2] = 0.  This has three real roots with just one
        satisifying q >= 0.
    '''
    a = l0 + l1
    b = l0 * l1
    c = -b * rs[2]  # z2
    # cubic equation z**3 + a*z**2 + b*z + c = 0
    b -= fdot(rs, l1, l0, a)
    a -= rt
    _r = b > 0
    if _r:
        a, b = b, a
        c  = _1_over(c)
        a *= c
        b *= c
    # see https://dlmf.nist.gov/1.11#iii
    p = (b * _3_0 - a**2) / _3_0
    t = -p / _3_0  # A / 4
    if t > 0:
        q = (a**3 * _2_0 - a * b * _9_0 + c * _27_0) / _27_0
        # switch to https://dlmf.nist.gov/4.43
        s = -q**2 - p**3 * _4_0 / _27_0
        p = sqrt(s) if s > 0 else _0_0
        s, c = sincos2(atan2(q, p) / _3_0)  # alp
        t = (c * _SQRT3 - s) * sqrt(t)
    else:
        t = _0_0
    t -= a / _3_0
    return _1_over(t) if _r else t


def _d_t(triax, llk):
    '''(INTERNAL) Helper.
    '''
    if llk in LLK._CENTRICS:
        d_t = True,  None
    elif llk in LLK._DETICS:
        d_t = True,  triax._a2b2c23
    elif llk in LLK._METRICS:
        d_t = False, triax._abc3
    else:
        raise TriaxialError(llk=llk)
    return d_t


def _h_llk_unit_name(height, h=None, llk=LLK.ELLIPSOIDAL, unit=Radians, **name):
    '''(INTERNAL) Helper, C{h} for backward compatibility.
    '''
    if llk is None:
        llk = LLK.ELLIPSOIDAL
    elif llk not in LLK._NOIDAL:  # or llk._X
        raise TriaxialError(llk=llk)
    if h is None:
        h = height
    return h, llk, unit, name


def _solve(_fp2, p, pscale, **n):
    '''(INTERNAL) Solve _fp2(p) = 0
    '''
    dt  = _N_1_0
    pt  = _EPSqrt * pscale
    _P2 =  Fsum(p).fsum2_
    for i in range(_MAXIT):
        fv, fp = _fp2(p, **n)
        if not (fv > _TOL2):
            break
        p, d = _P2(-fv / fp)  # d is positive
        if i and d <= dt and (fv <= EPS8 or
                 d <= (max(pt, p) * _TOL)):
            break
        dt = d
    else:
        t = Fmt.no_convergence(d, min(dt, pt))
        raise TriaxialError(_fp2.__name__, p, txt=t)
    return p


def _txtz2(cb, so, k, kp):
    '''(INTERNAL) Helper.
    '''
    return hypot(cb * k, kp), hypot(k, so * kp)


def _v_h_llk_name(x_ct, y=None, z=None, **h_llk_name):
    '''(INTERNAL) Helper.
    '''
    if y is z is None and isinstance(x_ct, Cartesian5Tuple):

        def _v_h_llk_name(h=x_ct.h, llk=x_ct.llk, **name):
            v = Vector3d(*x_ct.xyz3, **name)
            return v, h, llk, name
    else:
        def _v_h_llk_name(h=INT0, llk=None, **name):  # PYCHOK redef
            v = _otherV3d_(x_ct, y, z)
            return v, h, llk, name

    return _v_h_llk_name(**h_llk_name)


def _v_h_llk_name_NOIDAL(x_ct, y, z, **h_llk_name):
    '''(INTERNAL) Helper for methods C{reverseBetOmg} and C{-Alp}.
    '''
    v, h, llk, name = _v_h_llk_name(x_ct, y, z, **h_llk_name)
    if h or llk not in LLK._NOIDAL:  # or llk._X
        kwds = dict(x_ct=x_ct) if y is z is None else \
               dict(x=x_ct, y=y, z=z)
        raise TriaxialError(h=h, llk=llk, **kwds)
    return v, h, (LLK.ELLIPSOIDAL if llk is None else llk), name


class Triaxial3s(_TriaxialsBase):
    '''(INTERNAL) L{Triaxial3} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    _Triaxial = Triaxial3

Triaxial3s = Triaxial3s(Triaxial3, Triaxial3B)  # PYCHOK singleton
'''Some pre-defined L{Triaxial3}s, like L{Triaxials<triaxials.triaxial5.Triaxials>}.'''
Triaxial3s._assert()

if __name__ == _DMAIN_:
    # __doc__ of this file, force all into registry
    from pygeodesy.internals import _pregistry
    _pregistry(Triaxial3s)


# **) MIT License
#
# Copyright (C) 2025-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
