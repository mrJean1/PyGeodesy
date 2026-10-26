
# -*- coding: utf-8 -*-

u'''A pure Python implementation of 2 ellipsoidal and 3 spherical variants of I{Hotine}'s
U{Oblique Mercator<https://WikiPedia.org/wiki/Oblique_Mercator_projection>} projection from
U{Snyder Chapter 9<https://pubs.USGS.gov/pp/1395/report.pdf>} page 66ff, on-line U{Neacsu
<https://Neacsu.net/geodesy/snyder/3-cylindrical/sect_9>} and U{IOGP section 3.2.4
<https://www.IOGP.org/wp-content/uploads/2019/09/373-07-02.pdf>} page 59ff.

@see: U{Oblique Mercator Usage<https://Neacsu.net/geodesy/snyder/3-cylindrical/sect_9/#usage>}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _copysign,  _ALL_LAZY, typename
from pygeodesy.constants import NAN, PI_2, R_M, _0_0, _0_5, _1_0, _2_0, _90_0, \
                               _EPS2e4 as _TOL, _isfinite, isnear0, _1_over, _over
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _ValueError, _xkwds
from pygeodesy.fmath import fdot, fdot_
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import _C_, _COMMA_, _COMMASPACE_, _lat_, _lon_, _scale_, \
                              _x_, _y_, _0_, _1_, _2_
# from pygeodesy.lazily import _ALL_LAZY  # from .basics
from pygeodesy.named import _NamedBase, _NamedTuple, _namethod, _Pass
from pygeodesy.namedTuples import LatLon2Tuple, Vector2Tuple
from pygeodesy.props import Property_RO, property_doc_, property_RO
from pygeodesy.streprs import Fmt, pairs, unstr
from pygeodesy.units import Azimuth, Degrees, Int, Meter, Lamd, Lat, _LatLon2, Lon, \
                            Phid, Radius, Scalar
from pygeodesy.utily import asin1, atan1, atan2, sincos2, sincos2_, wrapPI as _PI

from math import cos, cosh, degrees, e as _base_e, fabs, log, radians, \
                 sin, sinh, sqrt, tan, tanh

__all__ = _ALL_LAZY.hom
__version__ = '26.10.24'


class _HotineObliqueBase(_NamedBase):
    '''(INTERNAL) Base class for all C{HotineOblique*} variants.
    '''
    _a0 = _b0 = _0_0  # origin (phi0, lam0)
    _ap = _bp = _0_0  # pole (phip, lamp)
    _args     = ()    # __init__ args
    _g0 = _zC =  NAN  # RSO angle, azimuth
    _k0       = _1_0  # central scale factor
    _lat0     = _0_0  # origin latitude (degrees) NAN
    _x0 = _y0 = _0_0  # false east-, northing

    def __init__(self, unused, k0=_1_0, lat0=None,
                               x0=_0_0, y0=_0_0, **name):
        if name:
            self.name = name
        if k0 != 1:
            self._k0 = Scalar(k0=k0)
#       self._Rk0 = self.R * self.k0
        if lat0 is not None:
            self._lat0 = Lat(lat0=lat0)
            self._a0   = Phid(self._lat0)
            sa0, ca0   = self._sca0
            if isnear0(sa0) or isnear0(ca0):
                raise self._Error(lat0=self._lat0)
        if x0:
            self._x0 = Meter(x0=x0)
        if y0:
            self._y0 = Meter(y0=y0)

    @Property_RO
    def _args_kwds(self):
        kwds = self._kwds()
        if self.name:
            kwds = _xkwds(kwds, name=self.name)
        return unstr(None, *self._args, **kwds)

    @Property_RO
    def azimuthC(self):
        '''Get the center line azimuth (C{degrees180}) or C{NAN}.
        '''
        z = self._zC
        if _isfinite(z):
            z = Azimuth(azimuth=degrees(z), _180=True)  # Snyder page 275, 𝛼c
        return z

    def _Error(self, **kwds):
        return HotineError(txt=str(self), **kwds)

    @Property_RO
    def gamma0(self):
        '''Get the Rectified Skew Orthomorphic (RSO) angle (C{degrees}) or C{NAN}.
        '''
        g = self._g0
        if _isfinite(g):
            g = Degrees(gamma0=degrees(g))
        return g

    @property_RO
    def k0(self):
        '''Get the scale factor at origin (C{scalar}).
        '''
        return self._k0

    def _kwds(self, **kwds):
        if self.k0 != 1:
            kwds.update(k0=self.k0)
        if self._lat0:
            kwds.update(lat0=self._lat0)
        if self._x0:
            kwds.update(x0=self._x0)
        if self._y0:
            kwds.update(y0=self._y0)
        return kwds

    @property_RO
    def latlon0(self):
        '''Get the origin, other pole (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(self._lat0, self._lon0, name=_0_)

    @Property_RO
    def _lon0(self):
        return degrees(self._b0)

    def _philam2(self, lat_latlon, lon):
        lat,  lon  = \
        self._args = _LatLon2(lat_latlon, Error=HotineError) if lon is None else (
                          Lat(lat_latlon), Lon(lon))
        return Phid(lat), Lamd(lon)

    def _philam4(self, latlon1, latlon2):
        a1, b1 = self._philam2(latlon1, None)
        t = self._args
        a2, b2 = self._philam2(latlon2, None)
        self._args = t, self._args
        return a1, b1, a2, b2

    def _philamalp3(self, latlon, lon, azimuth):
        a, b = self._philam2(latlon, lon)
        z    = Azimuth(azimuth)  # == % 360
        self._args += z,
        self._zC    = z = radians(z)
        self._b0    = b  # 𝜆C
        self._a0    = a  # 𝜑C
        return a, b, z

    @Property_RO
    def _sca0(self):  # 𝜑0
        return sincos2(self._a0)

#   @Property_RO
#   def _scb0(self):  # 𝜆0
#       return sincos2(self._b0)

    @Property_RO
    def _scap(self):  # 𝜑p
        return sincos2(self._ap)

#   @Property_RO
#   def _scbp(self):  # 𝜆p
#       return sincos2(self._bp)

    @Property_RO
    def _scg0(self):  # 𝛾0
        return sincos2(self._g0)

    @Property_RO
    def _sczC(self):  # 𝛼C or 𝛽C
        return sincos2(self._zC)

    def toStr(self, **kwds):
        '''Get this C{HotineOblique*} instance as C{str}.
        '''
        t = self._args_kwds
        return _COMMASPACE_(t, *pairs(kwds))

    @property_RO
    def xy0(self):
        '''Get the false east- and northing (L{Vector2Tuple}, C{meter}).
        '''
        return Vector2Tuple(self._x0, self._y0, name=_0_)


class _HotineOblique123(_HotineObliqueBase):
    '''(INTERNAL) Spherical C{HotineOblique1}, C{-2} and C{-3} base class.
    '''
    _Rk0    = \
    _radius = R_M

    def __init__(self, ap, bp, radius, **k0_lat0_x0_y0_name):
        self._ap = ap
        self._bp = bp = _PI(bp)
        self._b0 =      _PI(bp + PI_2)

        _HotineObliqueBase.__init__(self, None, **k0_lat0_x0_y0_name)

        if radius  != _HotineOblique123._radius or \
           self.k0 != _HotineObliqueBase._k0:
            self._radius = R = Radius(radius)
            self._Rk0    = R * self.k0

    def forward(self, lat_latlon, lon=None, **name):
        '''Project a geodetic lat- and longitude to Hotine east- and northing.

           @arg lat_latlon: Latitude of the point (scalar C{degrees90} or C{LatLon},
                    2-tuple C{(lat, lon)} or comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Longitude of the point (C{degrees180}), required if B{C{lat_latlon}}
                     is scalar, ignored otherwise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Hotine8Tuple}C{(x, y, u, v, lat, lon, scale, hom)} with falsed
                    east- and northing C{x} and C{y} and "unrectified" C{u} and C{v}
                    in C{meter}, geodetic C{lat}- and C{lon}gitude in degrees, C{scale}
                    factor (aka C{k}) and this C{HotineOblique*} Mercator instance C{hom}.
        '''
        lat, lon = _LatLon2(lat_latlon, lon, HotineError)

        d = Lamd(lon) - self._b0
        sa, ca, sd, cd = sincos2_(Phid(lat), d)
        sp, cp = self._scap  # phi-p

        A = fdot_(sp, sa, -cp, ca * sd)
        u = fdot_(cp, sa,  sp, ca * sd)  # / ca
        u = atan2( u, ca * cd)  * self._Rk0
        v = log(_d_s(_1_0, -A)) * self._Rk0 * _0_5

        k = _1_0 - A**2
        k = (self.k0 / sqrt(k)) if k > 0 else NAN
        x =  u + self._x0
        y =  v + self._y0
        return Hotine8Tuple(x, y, u, v, lat, lon, k, self,
                            name=_namethod(self.forward, **name))

    def _kwds(self, **unused):
        return _HotineObliqueBase._kwds(self, radius=self.radius)

    @Property_RO
    def latlonP(self):
        '''Get the main pole (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(degrees(self._ap), degrees(self._bp), name='P')

    @property_RO
    def radius(self):
        '''Get the the earth radius (C{meter}, conventionally).
        '''
        return self._radius

    def reverse(self, x, y, **name):
        '''Project an Hotine east- and northing to geodetic lat- and longitude.

           @arg x: Easting of the point (C{meter}).
           @arg y: Northing of the point (C{meter}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Hotine8Tuple}C{(x, y, u, v, lat, lon, scale, hom)} with geodetic
                    C{lat}- and C{lon}gitude in degrees, C{scale} factor (aka C{k}),
                    falsed east- and northing C{x} and C{y} and "unrectified" C(u) and
                    C{v} in C{meter} and this C{HotineOblique*} Mercator instance C{hom}.
        '''
        t = _1_over(self._Rk0)
        u = x - self._x0
        v = y - self._y0
        e = u * t
        n = v * t

        sx, cx     = sincos2(e)
        sy, cy, ty = sinh(n), cosh(n), tanh(n)
        sp, cp     = self._scap

        t = _over(sx, cy)
        a =  asin1(   fdot_(sp, ty,  cp, t))
        b = _PI(atan2(fdot_(sp, sx, -cp, sy), cx) + self._b0)
        return Hotine8Tuple(x, y, u, v, degrees(a), degrees(b), self.k0, self,
                            name=_namethod(self.reverse, **name))


class HotineOblique1(_HotineOblique123):
    '''Spherical variant 1: Hotine Oblique Mercator projection given
       C{two points}.
    '''

    def __init__(self, latlon1, latlon2, radius=R_M, **k0_lat0_x0_y0_name):
        '''New I{spherical} L{HotineOblique1} Mercator, variant 1.

           @arg latlon1: Center line point (C{LatLon}, C{2-tuple (lat, lon)} or
                         a comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg latlon2: Other center line point (C{LatLon}, C{2-tuple (lat, lon)}
                         or a comma separated I{"latDMS, lonDMS"} C{str}).
           @kwarg radius: Earth radius (C{meter}, conventionally).
           @kwarg k0_lat0_x0_y0_name: Central scale factor C{B{k0}=1.0} (C{scalar}), origin
                     latitude C{B{lat0}=0.0} (C{degrees}), false easting C{B{x0}=0.0} and
                     false northing C{B{y0}=0.0} (C{meter}) and C{B{name}=NN} (C{str}).
        '''
        a1, b1, a2, b2 = self._philam4(latlon1, latlon2)

        sa1, ca1, sb1, cb1 = sincos2_(a1, b1)
        sa2, ca2, sb2, cb2 = sincos2_(a2, b2)
        ca1sa2 = ca1 * sa2
        ca2sa1 = ca2 * sa1
        bp = atan2(fdot_(ca1sa2, cb1,  -ca2sa1, cb2),
                   fdot_(ca2sa1, sb2,  -ca1sa2, sb1))
        ap = atan1(-cos(bp - b1), tan(a1))
        _HotineOblique123.__init__(self, ap, bp, radius, **k0_lat0_x0_y0_name)

    @property_RO
    def latlon1(self):
        '''Get the 1st center line point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[0], name=_1_)

    @property_RO
    def latlon2(self):
        '''Get the 2nd center line point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[1], name=_2_)


class HotineOblique2(_HotineOblique123):
    '''Spherical variant 2: Hotine Oblique Mercator projection given the C{central point}
       and C{azimuth}.
    '''

    def __init__(self, lat_latlon, lon=None, azimuth=_0_0, radius=R_M, **k0_lat0_x0_y0_name):
        '''New spherical L{HotineOblique2} Mercator, variant 2.

           @arg lat_latlon: Center latitude (scalar C{degrees90} or any C{LatLon}, 2-tuple
                    C{(lat, lon)} or a comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Center longitude (C{degrees180}), required if B{C{lat_latlon}} is scalar,
                     =ignored otherwise.
           @kwarg azimuth: Azimuth of the center line (C{degrees360} from North).
           @kwarg radius: Earth radius (C{meter}, conventionally).
           @kwarg k0_lat0_x0_y0_name: Central scale factor C{B{k0}=1.0} (C{scalar}), origin
                     latitude C{B{lat0}=0.0} (C{degrees}), false easting C{B{x0}=0.0} and
                     false northing C{B{y0}=0.0} (C{meter}) and C{B{name}=NN} (C{str}).
        '''
        aC, bp, zC     = self._philamalp3(lat_latlon, lon, azimuth)
        sa, ca, sz, cz = sincos2_(aC, zC)  # == ._sca0, ._sczC

        ap  = asin1( sz * ca)
        bp += atan2(-cz, -sa * sz)  # PYLINT E1130
        _HotineOblique123.__init__(self, ap, bp, radius, **k0_lat0_x0_y0_name)

    @property_RO
    def latlonC(self):
        '''Get the center point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[:2], name=_C_)


class HotineOblique3(_HotineOblique123):
    '''Spherical variant 3: Hotine Oblique Mercator projection given the C{main pole}.
    '''

    def __init__(self, lat_latlon, lon=None, radius=R_M, **k0_lat0_x0_y0_name):
        '''New spherical L{HotineOblique3} Mercator, variant 3.

           @arg lat_latlon: Main pole latitude (scalar C{degrees90} or any C{LatLon},
                    2-tuple C{(lat, lon)} or a comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Main pole longitude (C{degrees180}), required if B{C{lat_latlon}} is
                     scalar, ignored otherwise.
           @kwarg radius: Earth radius (C{meter}, conventionally).
           @kwarg k0_lat0_x0_y0_name: Central scale factor C{B{k0}=1.0} (C{scalar}), origin
                     latitude C{B{lat0}=0.0} (C{degrees}), false easting C{B{x0}=0.0} and
                     false northing C{B{y0}=0.0} (C{meter}) and C{B{name}=NN} (C{str}).
        '''
        ap, bp = self._philam2(lat_latlon, lon)
        _HotineOblique123.__init__(self, ap, bp, radius, **k0_lat0_x0_y0_name)


class _HotineObliqueAB(_HotineObliqueBase):
    '''(INTERNAL) Ellipsoidal C{HotineObliqueA} and C{-B} base class.
    '''
    _ABEG_BuCz = ()  # 7-tuple
    _earth     = _WGS84
    _useRSO    =  False

    def __init__(self, earth, useRSO, **k0_lat0_x0_y0_name):

        if earth is not _HotineObliqueAB._earth:
            self._earth = _ellipsoidal_datum(earth)
        if useRSO:
            self.useRSO = useRSO

        _HotineObliqueBase.__init__(self, None, **k0_lat0_x0_y0_name)

    @Property_RO
    def _ABD21EFt0(self):  # 8-tuple
        sa0, ca0 = self._sca0
        E  =  self.ellipsoid
        e2 =  E.e2not
        r  = _1_0 - e2
        B  = _over(ca0**4 * e2, r)
        B  =  sqrt(_1_0 + B)
        s  = _1_0 - sa0**2 * e2
        r  =  sqrt(r) * B
        A  = _over(r * E.a * self.k0, s)
        a0 =  self._a0
        D  =  F = _over(r, sqrt(s) * ca0)
        D2 =  D**2
        if D2 > 1:
            D21 =  sqrt(D2 - _1_0)
            F  += _copysign(D21, a0)  # Snyder F, NOT!
        else:  # PYCHOK no cover
            D21 = _0_0
        t0 = _tesB(E.e, sa0, 1, a0)  # 9-13, not 9-13a!
        E  =  pow(t0, B) * F  # Snyder E, 9-15 & 9-36 IOGP H
        return A, B, D, D2, D21, E, F, t0

    @property_RO
    def ABmap(self):
        '''Get 14 map parameters C{(A, B, D, D2, E, G, uC, zC,
           a0, b0, g0, t0, x0, y0)} (L{HotineAB14Tuple}).
        '''
        A, B, D, D2, _, E, _, t0 = self._ABD21EFt0
        _, _, _, G, _, uC, _     = self._ABEG_BuCz  # PYLINT W0632
        return HotineAB14Tuple(A, B, D, D2, E, G, uC, self._zC,
                               self._a0, self._b0, self._g0, t0,
                               self._x0, self._y0, name=typename(self))

#   def _debug(self, where, locals__, **kwds):
#       d = dict(a0=self._a0, b0=self._b0, latlon0=self.latlon0,
#                g0=self._g0, k0=self.k0,  zC=self._zC)
#       d.update(locals__)
#       d.update(kwds)
#       n =  self.methodname(where)
#       s = _COMMA_ + _NL_ + _SPACE_ * (len(n) + 1)
#       print(n, s.join(pairs(d)))

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

    def forward(self, lat_latlon, lon=None, **name):
        '''Project a geodetic lat- and longitude to Hotine east- and northing.

           @arg lat_latlon: Latitude of the point (scalar C{degrees90} or C{LatLon},
                    2-tuple C{(lat, lon)} or comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Longitude of the point (C{degrees180}), required if B{C{lat_latlon}}
                     is scalar, ignored otherwise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Hotine8Tuple}C{(x, y, u, v, lat, lon, scale, hom)} with falsed east-
                    and northing C{x} and C{y} and "unrectified" C{u} and C{v} in C{meter},
                    geodetic C{lat}- and C{lon}gitude in degrees, C{scale} factor (aka
                    C{k}) and this C{HotineOblique*} Mercator instance C{hom}.
        '''
        lat, lon = _LatLon2(lat_latlon, lon, HotineError)

        A, B, Q, _, A_B, uC, zCE = self._ABEG_BuCz  # PYLINT W0632

        a  =  Phid(lat)
        b  =  Lamd(lon)
        dB = _PI(b - self._b0) * B
        E  =  self.ellipsoid
        sg0,  cg0 = self._scg0

        sa, ca, V, cdB = sincos2_(a, dB)  # phi
        if ca:  # PYLINT W0125
            Q   *= _tesB(E.e, sa, -B, a)  # PYLINT E1120 Snyder 9-25 ...
            S, T = _ST_Q(Q)
            U    = _over(fdot_(S, sg0, -V, cg0), T)  # PYLINT E1130

            v = _d_s(_1_0, U)
            if _isfinite(v) and v > 0:  # U +/-1
                v = log(v) * A_B * _0_5  # ... Snyder 9-30
#           elif ca:  # PYLINT W0125
#               pass  # 0.0 or INF
#           else:  # lat +/-90
#               v = _copysign(self._g0 * _0_5, -sa)  # PYLINT E1130
#               v =  log(tan(PI_4 + v)) * A_B

            if cdB:  # PYLINT W0125
                u = fdot_(S, cg0,   V, sg0)
                u = atan2(u, cdB) * A_B  # Snyder 9-31 XXX atan1
            else:  # PYCHOK no cover
                u = A * dB  # Note after 9-31

            if uC or zCE:  # variant B, IOGP page 62
                u -= (uC if b < self._b0 else  # ABS(uC) * SIGN(a0) * SIGN(b0 - b)
                    (-uC if b > self._b0 else u)) if zCE else uC
        else:  # PYCHOK no cover
            u =  a * A_B  # lat +/-90 idem
            v = _0_0

        x = _1_0 - sa**2 * E.e2not
        if x > 0:
            k = E.a * ca * cdB  # Snyder 9-32
            k = Scalar(k=cos(u / A_B) * sqrt(x) * A / k) if k and A_B else NAN
        else:  # PYCHOK no cover
            k = NAN if x < 0 else _0_0
        x, y = self._rectify2(u, v, self._x0, self._y0)

#       self._debug(self.forward, locals(), t=_tesB(E.e, sa, 1, a))
        return Hotine8Tuple(x, y, u, v, lat, lon, k, self,
                            name=_namethod(self.forward, **name))

    def _kwds(self, **unused):
        return _HotineObliqueBase._kwds(self, earth=self.ellipsoid, useRSO=self.useRSO)

    def _rectify2(self, e, n, uC, v0):
        s, c = self._scg0 if self.useRSO else self._sczC  # PYLINT W0125
        return (fdot_(n, c,  e, s, start=uC),  # 9-40
                fdot_(e, c, -n, s, start=v0))  # 9-41

    def reverse(self, x, y, triPhi=0, **name):
        '''Project an Hotine east- and northing to geodetic lat- and longitude.

           @arg x: Easting of the point (C{meter}).
           @arg y: Northing of the point (C{meter}).
           @kwarg triPhi: Number of C{phi} iterations (C{int}), set C{B{triPhi}=0}
                          to use the Horner-like C{e2-sin(X2)} series.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Hotine8Tuple}C{(x, y, u, v, lat, lon, scale, hom)} with geodetic
                    C{lat}- and C{lon}gitude in degrees, C{scale} factor (aka C{k}),
                    falsed east- and northing C{x} and C{y} and "unrectified" C{u} and
                    C{v} in C{meter} and this C{HotineOblique*} Mercator instance C{hom}.
        '''
        _, B, E, _, A_B, uC, _ = self._ABEG_BuCz  # PYLINT W0632

        u, v     = self._rectify2(x - self._x0, y - self._y0, uC, _0_0)
        V, cuB_A = sincos2(_over(u, A_B))
        sg0, cg0 = self._scg0

        Q    =  pow(_base_e, _over(-v, A_B))
        S, T = _ST_Q(Q)
        U    = _over(fdot_(V, cg0,  S, sg0), T)
        t, i = _d_s(_1_0, -U), None
        if _isfinite(t) and t > 0:
            t =  pow(E / sqrt(t), _1_over(B))
            X = _phit(t)  # Snyder chi, 7-13
            n =  Int(triPhi=triPhi)
            e =  self.ellipsoid.e
            if e and n > 0:
                e_2 = e * _0_5
                for i in range(n):  # 4..5
                    a = _d_s(_1_0, sin(X) * e)
                    a = _phit(pow(a, e_2) * t)  # 7-9
                    d =  fabs(a - X)
                    if d < _TOL:
                        i += 1
                        break
                    X = a
                else:  # PYCHOK no cover
                    t = unstr(self.reverse, x, y, triPhi=n)
                    e = Fmt.no_convergence(d, _TOL)
                    raise HotineError(e, txt=t)

            else:  # Horner-like e2-sin(X2) series, 3-5
                sX4 = [sin(X * _) for _ in range(2, 9, 2)]
                a   = fdot(self._Xe5, X, *sX4)

            b =  fdot_(S, cg0, -V, sg0)
            b =  atan2(b, cuB_A)
            b = _PI(self._b0 - _over(b, B))

            lat, lon =  degrees(a), degrees(b)
        else:  # PYCHOK no cover  a = +/-PI_2, b = self._b0
            lat, lon = _copysign(_90_0, U), self._lon0

#       self._debug(self.reverse, locals())
        return Hotine8Tuple(x, y, u, v, lat, lon, self.k0, self, iteration=i,
                            name=_namethod(self.reverse, **name))

    @property_doc_(' apply Rectified Skew Orthomorphic (RSO) angle C{gamma0}, otherwise Oblique Mercator C{azimuthC} (C{bool}).')
    def useRSO(self):
        return self._useRSO

    @useRSO.setter  # PYCHOK setter!
    def useRSO(self, RSO):
        self._useRSO = bool(RSO)

    @Property_RO
    def _Xe5(self):  # Snyder page 75, 3-5
        e2 = self.ellipsoid.e2not
        e4 = e2 * e2
        e6 = e4 * e2
        return (_1_0, e2 * (_0_5 + e2 * (5 / 24 + e2 * ( 1 / 12  + e2 * (  13 / 360)))),
                      e4 * (             7 / 48 + e2 * (29 / 240 + e2 * ( 811 / 11520))),
                      e6 * (                             7 / 120 + e2 * (  81 / 1120)),
                      e6 *                                         e2 * (4279 / 161280))


class HotineObliqueA(_HotineObliqueAB):
    '''Ellipsoidal variant A: Hotine Oblique Mercator projection given C{two points}.
    '''

    def __init__(self, latlon1, latlon2, earth=_WGS84, useRSO=False, **k0_lat0_x0_y0_name):
        '''New ellipsoidal L{HotineObliqueA} Mercator, variant A.

           @arg latlon1: Center line point (C{LatLon}, 2-tuple C{(lat, lon)} or
                         comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg latlon2: Other center line point (C{LatLon}, 2-tuple C{(lat, lon)} or
                         comma-separated I{"latDMS, lonDMS"} C{str}).
           @kwarg earth: Earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg useRSO: If C{True} apply Rectified Skew Orthomorphic (RSO) angle C{gamma0},
                          otherwise Oblique Mercator C{azimuthC} (C{bool}).
           @kwarg k0_lat0_x0_y0_name: Central scale factor C{B{k0}=1.0} (C{scalar}), origin
                     latitude C{B{lat0}=0.0} (C{degrees}), false easting C{B{x0}=0.0} and
                     false northing C{B{y0}=0.0} (C{meter}) and C{B{name}=NN} (C{str}).
        '''
        a1, b1, a2, b2 = self._philam4(latlon1, latlon2)
        _HotineObliqueAB.__init__(self, earth, useRSO, **k0_lat0_x0_y0_name)

        (lat1, _), (lat2, _) = self._args
        if isnear0(a1 - a2):  # phi1 == phi2
            raise self._Error(lat1=lat1, lat2=lat2)
        sa1, ca1, sa2, ca2 = sincos2_(a1, a2)
        if isnear0(sa1) or isnear0(ca1):  # phi1 +/-PI_2 or 0
            raise self._Error(lat1=lat1)
        elif isnear0(ca2 + PI_2):  # phi2 == -PI_2
            raise self._Error(lat2=lat2)

        A, B, D, _, _, E, _, _ = self._ABD21EFt0
        e  =  self.ellipsoid.e
        H  = _tesB(e, sa1, B, a1)  # Snyder H
        L  = _tesB(e, sa2, B, a2)
        G  = (_over(E, H) - _over(H, E)) * _0_5
        J  = _d_s(E**2, L * H)
        P  = _d_s(L, H)

        b0 =  atan1(tan(B * _PI(b1 - b2) * _0_5) * J, P)  # atan2?
        b0 = (b1 + b2) * _0_5 - _over(b0, B)
        g0 =  atan1(sin(B * _PI(b1 - b0)), G)  # atan2?
        zC =  asin1(sin(g0) * D)
        uC = _0_0

        self._b0, self._g0, self._zC = b0, g0, zC
        self._ABEG_BuCz = A, B, E, G, _over(A, B), uC, False
#       self._debug(HotineObliqueA.__init__, locals(), gamma0=self.gamma0, zCd=self.azimuth)

    @property_RO
    def latlon1(self):
        '''Get the 1st center line point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[0], name=_1_)

    @property_RO
    def latlon2(self):
        '''Get the 2nd center line point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[1], name=_2_)


class HotineObliqueB(_HotineObliqueAB):
    '''Ellipsoidal variant B: Hotine Oblique Mercator projection given the C{central point}
       and C{azimuth}.
    '''

    def __init__(self, lat_latlon, lon=None, azimuth=_0_0, earth=_WGS84, useRSO=False, **k0_lat0_x0_y0_name):
        '''New ellipsoidal L{HotineObliqueB} Mercator, variant B.

           @arg lat_latlon: Center latitude (scalar C{degrees90} or C{LatLon}, 2-tuple
                    C{(lat, lon)} or comma-separated I{"latDMS, lonDMS"} C{str}).
           @arg lon: Center longitude (C{degrees180}), required if B{C{lat_latlon}} is scalar,
                     ignored otherwise.
           @kwarg azimuth: Azimuth of the center line (C{degrees360} from North).
           @kwarg earth: Earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg useRSO: If C{True} apply Rectified Skew Orthomorphic (RSO) angle C{gamma0},
                          otherwise Oblique Mercator C{azimuthC} (C{bool}).
           @kwarg k0_lat0_x0_y0_name: Central scale factor C{B{k0}=1.0} (C{scalar}), origin
                     latitude C{B{lat0}=0.0} (C{degrees}), false easting C{B{x0}=0.0} and
                     false northing C{B{y0}=0.0} (C{meter}) and C{B{name}=NN} (C{str}).
        '''
        aC, bC, zC = self._philamalp3(lat_latlon, lon, azimuth)
        _HotineObliqueAB.__init__(self, earth, useRSO, **k0_lat0_x0_y0_name)
        szC, czC = self._sczC

        A, B, D, _, D21, E, F, _ = self._ABD21EFt0
        A_B = _over(A, B)
        G   = (F - _1_over(F)) * _0_5

        g0  =  asin1(_over(szC, D))  # == IOGP
        b0  =  bC - _over(asin1(tan(g0) * G), B)  # == IOGP

        zCE = isnear0(zC - PI_2)  # IOGP page 61, zC == 90?
        uC  = (A * (bC - b0)) if zCE else (
               atan2(D21, czC) * A_B)
        uC  = _copysign(uC, aC)  # == ABS(uC) * SIGN(𝜑C)

        self._b0, self._g0 = _PI(b0), g0

        self._ABEG_BuCz = A, B, E, G, A_B, uC, zCE
#       self._debug(HotineObliqueB.__init__, locals(), gamma0=self.gamma0)

    @property_RO
    def latlonC(self):
        '''Get the center point (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(*self._args[:2], name=_C_)


class HotineError(_ValueError):
    '''Error raised for C{HotineOblique*} issues.
    '''
    pass  # PYLINT W0107


class Hotine8Tuple(_NamedTuple):
    '''8-Tuple C{(x, y, u, v, lat, lon, scale, hom)} with (falsed) east- and northing C{x}
       and C{y} and "unrectified" C{u} and C{v} in C{meter}, geodetic C{lat}- and C{lon}gitude
       in degrees, C{scale} factor (aka C{k}) and C{HotineOblique*} Mercator instance C{hom}.

       @note: Except C{hom}, any attibute may be C{NAN} or non-finite for singulaties.
    '''
    _Names_ = (_x_,   _y_,   'u',   'v',   _lat_, _lon_, _scale_, 'hom')
    _Units_ = ( Meter, Meter, Meter, Meter, Lat,   Lon,  _Pass,   _Pass)


class HotineAB14Tuple(_NamedTuple):
    '''14-Tuple C{(A, B, D, D2, E, G, uC, zC, a0, b0, g0, t0, x0, y0)} of ellipsoidal map
       parameters I{only} for variants L{HotineObliqueA} and L{HotineObliqueB} with
       angles C{zC}, C{a0}, C{b0} and C{g0} in C{radians}.
    '''
    _Names_ = tuple('A,B,D,D2,E,G,uC,zC,a0,b0,g0,t0,x0,y0'.split(_COMMA_))
    _Units_ = (Scalar,) * 7 + (_Pass,) * 4 + (Scalar,) + (Meter,) * 2


def _d_s(a, b):  # diff / sum
    return _over(a - b, a + b)


def _phit(t):  # Snyder page 75, X 7-13
    return PI_2 - atan1(t) * _2_0  # 7-13


def _ST_Q(Q):
    _Q = _1_over(Q)
    S  = (Q - _Q) * _0_5
    T  = (Q + _Q) * _0_5
    return S, T


def _tesB(e, s, B, *a):  # t(e, s) power B!
    if a and B:  # ϕ, Snyder page 71, 9-13
        t  = _d_s(_1_0, e *  s)
        t  =  pow( t,  -e * _0_5)
        t *=  tan(_0_5 * (PI_2 - a[0]))
        t  =  pow( t, B)
#   elif s and B:  # Snyder page 71, 9-13a
#       t  = _d_s(_1_0, -s * e)
#       t  =  pow(t, e)
#       t *= _d_s(_1_0,  s)
#       t  =  pow(t, B * _0_5)  # == pow(sqrt(t), B)
    else:  # PYCHOK no cover
        t  = _1_0
    return t


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
