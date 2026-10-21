
# -*- coding: utf-8 -*-

u'''A pure Python implementation of I{Laborde}'s U{Oblique Mercator
<https://WikiPedia.org/wiki/Oblique_Mercator_projection>} projection U{IOGP section
3.2.4.2<https://www.IOGP.org/wp-content/uploads/2019/09/373-07-02.pdf>} page 64ff.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _copysign,  _ALL_LAZY
from pygeodesy.constants import PI_2, _0_0, _0_5, _1_0, _2_0, _3_0, _over, \
                               _EPS2e4 as _TOL
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _ValueError
from pygeodesy.fmath import fdot_, hypot
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import _C_, _COMMASPACE_, _lat_, _lon_, _scale_, \
                              _V_, _x_, _y_, _0_
# from pygeodesy.lazily import _ALL_LAZY  # from .basics
from pygeodesy.named import _NamedBase, _NamedTuple, _namethod, _Pass
from pygeodesy.namedTuples import LatLon2Tuple, PhiLam2Tuple, Vector2Tuple
from pygeodesy.props import Property_RO, property_RO, property_ROver
from pygeodesy.streprs import Fmt, pairs, unstr
from pygeodesy.units import Azimuth, Meter, Lamd, Lat, _LatLon2, Lon, \
                            Phid, Scalar
from pygeodesy.utily import asin1, atan1, sincos2, sincos2_, wrapPI as _PI

from math import atan, degrees, e as _base_e, fabs, log, radians, \
                 sin, sqrt, tan

__all__ = _ALL_LAZY.lom
__version__ = '26.10.20'

_TRIPS = 17  # per iteration loop


class LabordeOblique(_NamedBase):
    '''Ellipsoidal Laborde Oblique Mercator projection given the C{central point}
       and C{azimuth}.
    '''
    _aC    = _0_0  # phiC, radians
    _args  = ()  # 3-tuple
    _aS    = _0_0  # phiS, radians
    _bC    = _0_0  # lamC, radians
    _BCGR  = ()    # 4-tuple
    _d     = _0_0  # iteration delta
    _earth = _WGS84
    _kC    = _1_0  # central scale
    _x0    = _0_0  # false easting
    _y0    = _0_0  # false northing
    _zC    = _0_0  # central azimuth, radians

    def __init__(self, lat_latlon, lon=None, azimuth=_0_0, earth=_WGS84, kC=_1_0, x0=_0_0, y0=_0_0, **name):
        '''New ellipsoidal L{LabordeOblique} Mercator.

           @arg lat_latlon: Central latitude (scalar C{degrees90} or C{LatLon}, 2-tuple
                    C{(lat, lon)} or comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Central longitude (C{degrees180}), required if B{C{lat_latlon}} is scalar,
                     ignored otherwise.
           @kwarg azimuth: Azimuth of the center line (C{degrees360} from North).
           @kwarg earth: Earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg kC: Central scale factor (C{scalar}).
           @kwarg x0: False easting (C{meter}).
           @kwarg y0: False northing (C{meter}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @see: U{Roundtrip error<ttps://docs.Carmenta.com/pages/labordeobliquemercatorprojection.html>}
        '''
        self._args = lat, lon = _LatLon2(lat_latlon, lon, LabordeError)
        self._aC, self._bC    =  Phid(lat), Lamd(lon)

        zC = Azimuth(azimuth)
        self._args += zC,
        if zC:
            self._zC  = zC = radians(zC)
            szC, czC, s, c = sincos2_(zC, zC * _2_0)
            G = complex(_1_0 - c, s) / 12.0
            del s, c
        else:  # PYCHOK no cover
            czC =          _1_0
            szC = G = zC = _0_0
        if earth is not LabordeOblique._earth:
            self._earth = _ellipsoidal_datum(earth)

        if name:
            self.name = name
        if kC != 1:
            self._kC = Scalar(kC=kC)
        if x0:
            self._x0 = Meter(x0=x0)
        if y0:
            self._y0 = Meter(y0=y0)

        E =  self.ellipsoid
        r = _1_0 - E.e2not

        B = _over(czC**4 * E.e2not, r)
        B =  sqrt(B + _1_0)

        aS = asin1(_over(szC, B))
        C  = log(_tan_2(aS)) - _log_2(zC, E.e) * B

        R = _1_0 - szC**2 * E.e2not
        R = _over(sqrt(r) * E.a * kC, R)

        self._aS   = aS
        self._BCGR = B, C, G, R

    @Property_RO
    def _args_kwds(self):
        kwds = dict(earth=self.ellipsoid)
        if self.kC != 1:
            kwds.update(kC=self.kC)
        if self._x0:
            kwds.update(x0=self._x0)
        if self._y0:
            kwds.update(y0=self._y0)
        if self.name:
            kwds.update(name=self.name)
        return unstr(None, *self._args, **kwds)

    @Property_RO
    def azimuthC(self):
        '''Get the center line azimuth (C{degrees180}) or C{NAN}.
        '''
        return Azimuth(azimuth=degrees(self._zC), _180=True)

    @property_RO
    def earth(self):
        '''Get the earth model (L{Datum}).
        '''
        return self._earth

    @Property_RO
    def ellipsoid(self):
        '''Get the earth' ellipsoid (L{Ellipsoid}).
        '''
        return self.earth.ellipsoid

    def forward(self, lat_latlon, lon, **name):
        '''Project a geodetic lat- and longitude in C{degrees} to Laborde east- and northing.

           @arg lat_latlon: Latitude of the point (scalar C{degrees90} or C{LatLon},
                    2-tuple C{(lat, lon)} or comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Longitude of the point (C{degrees180}), required if B{C{lat_latlon}}
                     is scalar, ignored otherwise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Laborde8Tuple}C{(x, y, U, V, lat, lon, scale, lom)} with falsed
                    east- and northing C{x} and C{y} in C{meter}, "unrectified" C{U} and
                    C{V}, geodetic C{lat}- and C{lon}gitude in degrees, C{scale} factor
                    (aka C{k}) and this C{LabordeOblique*} Mercator instance C{lom}.
        '''
        lat, lon = _LatLon2(lat_latlon, lon, LabordeError)
        a,   b   =  Phid(lat), Lamd(lon)

        B, C, G, R = self._BCGR
        e = self.ellipsoid.e

        L = (b - self._bC) * B
        q = _log_2(a, e) * B + C
        _, P = _phiq2(q)
        L, P, U, V = _LPUV4(L, P, True, *self._scaS)

        H  = complex(-L, P)
        H += pow(H, 3) * G
        x  = H.imag * R + self._x0
        y  = H.real * R + self._y0
        return Laborde8Tuple(x, y, U, V, lat, lon, self.kC, self,
                             name=_namethod(self.forward, **name))

    @property_RO
    def kC(self):
        '''Get the scale factor at center (C{scalar}).
        '''
        return self._kC

    @property_RO
    def latlonC(self):
        '''Get the center point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[:2], name=_C_)

    @property_ROver
    def latlonParis(self):
        '''Get Paris Observatory lat- and Paris Meridian longitude (L{LatLon2Tuple}).
        '''
        # <https://Wikipedia.org/wiki/Paris_Observatory> <https://WikiPedia.org/wiki/Paris_meridian>
        return LatLon2Tuple(Lat('48 48 18.36'), Lon('2 20 14.025'), name='Paris')

    @property_RO
    def philamC(self):
        '''Get the center point (L{PhiLam2Tuple}).
        '''
        return _philam2(self.latlonC)

    @property_ROver
    def philamParis(self):
        '''Get Paris Observatory lat- and Paris Meridian longitude (L{PhiLam2Tuple}).
        '''
        return _philam2(self.latlonParis)

    def reverse(self, x, y, **name):
        '''Project an Laborde east- and northing to geodetic lat- and longitude.

           @arg x: Easting of the point (C{meter}).
           @arg y: Northing of the point (C{meter}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Laborde8Tuple}C{(x, y, U, V, lat, lon, scale, lom)} with geodetic
                    C{lat}- and C{lon}gitude in degrees, "unrectified" C{U} and C{V},
                    C{scale} factor (aka C{k}), falsed east- and northing C{x} and C{y}
                    and this C{LabordeOblique*} Mercator instance C{lom}.
        '''
        B, C, G, R = self._BCGR

        H0 =  complex(y - self._y0,
                      x - self._x0) / R
        P  =  pow( H0, 3) * G
        H_ = _over(H0, P + H0)
        for i in self._TRIPS(0, x, y):
            L =  pow(H_, 2) * G
            P =  pow(H_, 3) * G
            H = _over(P * _2_0 +  H0,
                      L * _3_0 + _1_0)
            P += H
            P -= H0
            if self._TOL(P.real):
                break
            H_ = H
        P = pow(_base_e, H.imag)
        P = atan(P) * _2_0 - PI_2
        L, P, U, V = _LPUV4(-H.real, P, False, *self._scaS)

        b     = _PI(_over(L, B) + self._bC)
        e     =  self.ellipsoid.e
        q, a_ = _phiq2(_over(P - C, B))
        for i in self._TRIPS(i, x, y):  # PYLINT W0631
            a = _1_se(-sin(a_), e) * q
            a = _phiq(a)
            if self._TOL(a - a_):
                break
            a_ = a

        return Laborde8Tuple(x, y, U, V, degrees(a), degrees(b), self.kC, self,
                             iteration=i, name=_namethod(self.reverse, **name))

    @Property_RO
    def _scaS(self):
        return sincos2(self._aS)

#   @Property_RO
#   def _sczC(self):
#       return sincos2(self._zC)

    def toStr(self, **kwds):
        '''Get this C{HotineOblique*} instance as C{str}.
        '''
        t = self._args_kwds
        return _COMMASPACE_(t, *pairs(kwds))

    def _TOL(self, d):
        self._d = d = fabs(d)
        return d < _TOL

    def _TRIPS(self, n, *xy):
        for i in range(n + 1, n + _TRIPS):
            yield i
        else:  # PYCHOK no cover
            t = unstr(self.reverse, *xy)
            e = Fmt.no_convergence(self._d, _TOL)
            raise LabordeError(e, txt=t)

    @property_RO
    def xy0(self):
        '''Get the false east- and northing (L{Vector2Tuple}, C{meter}).
        '''
        return Vector2Tuple(self._x0, self._y0, name=_0_)


class LabordeError(_ValueError):
    '''Error raised for C{LabordeOblique*} issues.
    '''
    pass  # PYLINT W0107


class Laborde8Tuple(_NamedTuple):
    '''8-Tuple C{(x, y, U, V, lat, lon, scale, lom)} with (falsed) east- and northing
       C{x} and C{y} in C{meter}, "unrectified" C{U} and C{V}, geodetic C{lat}- and
       C{lon}gitude in degrees, C{scale} factor (aka C{k}) and C{LabordeOblique*}
       Mercator instance C{lom}.

       @note: Except C{hom}, any attribute may be C{NAN} or non-finite for singulaties.
    '''
    _Names_ = (_x_,   _y_,   'U',    _V_,    _lat_, _lon_, _scale_, 'lom')
    _Units_ = ( Meter, Meter, Scalar, Scalar, Lat,   Lon,   Scalar, _Pass)


def _log_2(a, e):
    r  = _1_se(sin(a), e)
    r *= _tan_2(a)
    return log(r)


def _LPUV4(L, P, forward, saS, caS):
    sL, cL, sP, cP = sincos2_(L, P)
    V =                     cP * cL
    W = sP if forward else (cP * sL)
    U = fdot_(V, caS,  W, saS)
    V = fdot_(V, saS, -W, caS)
    if forward:
        W = cP * sL
    else:
        W = V
        V = sP
    d = hypot(U, V)
    if d:
        L = atan1(V, U + d) * _2_0
        P = atan1(W, d)
    else:  # PYCHOK no cover
        L = _0_0
        P = _copysign(PI_2, W)
    return L, log(_tan_2(P)), U, V


def _philam2(t):
    return PhiLam2Tuple(radians(t.lat), radians(t.lon), name=t.name)


def _phiq(q):
    return atan(q) * _2_0 - PI_2


def _phiq2(q):
    q = pow(_base_e, q)
    return q, _phiq(q)


def _1_se(s, e):  # diff / sum
    s *=  e
    r  = _over(_1_0 - s, _1_0 + s)
    return pow(r, e * _0_5)


def _tan_2(a):
    # tan(x/2) == sin(x) / (cos(x) + 1)
    #          == (1 - cos(x)) / sin(x)
    return tan((PI_2 + a) * _0_5)

# **) MIT License
#
# Copyright (C) 2026-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
