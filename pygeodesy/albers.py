
# -*- coding: utf-8 -*-

u'''Albers Equal-Area projections.

Classes L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth}, L{AlbersEqualAreaSouth} and
L{AlbersError}, a transcription of I{Charles Karney}'s C++ class U{AlbersEqualArea
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1AlbersEqualArea.html>}.

See also I{Albers Equal-Area Conic Projection} in U{John P. Snyder, "Map Projections
-- A Working Manual", 1987<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 98-106
and the Albers Conical Equal-Area examples on pp 291-294.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import copysign0, neg
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _ValueError, _xkwds
from pygeodesy.fmath import Fsum, fsum_, hypot, hypot1, sqrt3
from pygeodesy.interns import EPS0, NN, _EPSqrt as _TOL, _datum_, \
                             _EPS__2, _EPS__4, _gamma_, _lat_, \
                             _lat1_, _lat2_, _lon_, _no_, _scale_, \
                             _x_, _y_, _0_0, _0_5, _1_0, _2_0, \
                             _3_0, _90_0
from pygeodesy.karney import _diff182, _norm180
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass
from pygeodesy.props import deprecated_Property_RO, Property_RO, \
                            property_RO
from pygeodesy.streprs import Fmt
from pygeodesy.units import Bearing, Float_, Lat, Lat_, Lon, Lon_, \
                            Meter, Scalar_
from pygeodesy.utily import atand, atan2d, degrees360, sincos2, sincos2d

from math import atan, atan2, atanh, degrees, radians, sqrt

__all__ = _ALL_LAZY.albers
__version__ = '21.06.10'

_NUMIT  =  8  # XXX 4?
_NUMIT0 = 41  # XXX 21?
_TERMS  = 21  # XXX 9?
_TOL0   = sqrt3(_TOL)


class AlbersError(_ValueError):
    '''An L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
       L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth},
       L{AlbersEqualAreaSouth} or L{Albers7Tuple} issue.
    '''
    pass


def _Dsn(x, y, sx, sy):
    '''(INTERNAL) Divided differences, defined as M{Df(x, y) = (f(x) - f(y)) / (x - y)}
       with M{sn(x) = x / sqrt(1 + x^2)}: M{Dsn(x, y) = (x + y) / ((sn(x) + sn(y)) *
       (1 + x^2) * (1 + y^2))}.

       @see: U{W. M. Kahan and R. J. Fateman, "Sympbolic Computation of Divided
             Differences"<https://People.EECS.Berkeley.EDU/~fateman/papers/divdiff.pdf>},
             U{ACM SIGSAM Bulletin 33(2), 7-28 (1999)<https://DOI.org/10.1145/334714.334716>}
             and U{AlbersEqualArea.hpp
             <https://GeographicLib.sourceforge.io/html/AlbersEqualArea_8hpp_source.html>}.
    '''
    # sx = x / hypot1(x)
    d, t = _1_0, x * y
    if t > 0:
        s = sx + sy
        if s:
            d = (x + y) * ((sx * sy) / t)**2 / s
    else:
        t = x - y
        if t:
            d = (sx - sy) / t
    return d


def _Ks(**name_k):
    '''(INTERNAL) Scale C{B{k} >= EPS0}.
    '''
    return Scalar_(Error=AlbersError, low=EPS0, **name_k)  # > 0


def _Lat(**name_lat):
    '''(INTERNAL) Latitude C{-90 <= B{lat} <= 90}.
    '''
    return Lat_(Error=AlbersError, **name_lat)


def _Lon(**name_lon):
    '''(INTERNAL) Longitude C{-180 <= B{lon} <= 180}.
    '''
    return Lon_(Error=AlbersError, **name_lon)


def _tol(tol, x):
    '''(INTERNAL) Converge tolerance.
    '''
    return tol * max(_1_0, abs(x))


class _AlbersBase(_NamedBase):
    '''(INTERNAL) Base class for C{AlbersEqualArea...} projections.

       @see: I{Karney}'s C++ class U{AlbersEqualArea<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1AlbersEqualArea.html>}, method C{Init}.
    '''
    _datum     = _WGS84
    _iteration =  None
    _k0        =  None  # scale
    _k0n0      =  None  # (INTERNAL) k0 * no
    _k02       =  None  # (INTERNAL) k0**2
    _k02n0     =  None  # (INTERNAL) k02 * n0
    _lat0      =  None  # lat origin
    _lat1      =  None  # let 1st parallel
    _lat2      =  None  # lat 2nd parallel
    _m0        = _0_0   # if polar else sqrt(m02)
    _m02       =  None  # (INTERNAL) cached
    _n0        =  None  # (INTERNAL) cached
    _nrho0     = _0_0   # if polar else m0 * E.a
    _polar     =  False
    _qx        =  None  # (INTERNAL) cached
    _qZ        =  None  # (INTERNAL) cached
    _qZa2      =  None  # (INTERNAL) qZ * E.a**2
    _scxi0     =  None  # (INTERNAL) sec(xi)
    _sign      = +1
    _sxi0      =  None  # (INTERNAL) sin(xi)
    _txi0      =  None  # (INTERNAL) tan(xi)

    def __init__(self, sa1, ca1, sa2, ca2, k, datum, name):
        '''(INTERNAL) New C{AlbersEqualArea...} instance.
        '''
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)
        if name:
            self.name = name

        E = self.datum.ellipsoid
        b_a = E.b_a  # fm  = 1 - E.f
        e2  = E.e2
        e12 = E.e12  # e2m = 1 - E.e2

        self._qZ   = qZ = _1_0 + e12 * self._atanhee(1)
        self._qZa2 = qZ * E.a2
        self._qx   = qZ / (_2_0 * e12)

        c = min(ca1, ca2)
        if c < 0:
            raise AlbersError(clat1=ca1, clat2=ca2)
        polar = c < _EPS__2  # == 0
        # determine hemisphere of tangent latitude
        if sa1 < 0:  # and sa2 < 0:
            self._sign = -1
            # internally, tangent latitude positive
            sa1, sa2 = -sa1, neg(sa2)
        if sa1 > sa2:  # make phi1 < phi2
            sa1, sa2 = sa2, sa1
            ca1, ca2 = ca2, ca1
        if sa1 < 0:  # or sa2 < 0:
            raise AlbersError(slat1=sa1, slat2=sa2)
        # avoid singularities at poles
        ca1, ca2 = max(_EPS__2, ca1), max(_EPS__2, ca2)
        ta1, ta2 = sa1 / ca1, sa2 / ca2

        par1 = abs(ta1 - ta2) < _EPS__4  # ta1 == ta2
        if par1 or polar:
            C, ta0 = _1_0, ta2
        else:
            s1_qZ, C = self._s1_qZ_C2(ca1, sa1, ta1, ca2, sa2, ta2)

            ta0 = (ta2 + ta1) * _0_5
            Ta0 =  Fsum(ta0)
            tol = _tol(_TOL0, ta0)
            for self._iteration in range(1, _NUMIT0):
                ta02  = ta0**2
                sca02 = ta02 + _1_0
                sca0  = sqrt(sca02)
                sa0   = ta0 / sca0
                sa01  = sa0 + _1_0
                sa02  = sa0**2
                # sa0m = 1 - sa0 = 1 / (sec(a0) * (tan(a0) + sec(a0)))
                sa0m  = _1_0 / (sca0 * (ta0 + sca0))  # scb0^2 * sa0
                g  = (_1_0 + (b_a * ta0)**2) * sa0
                dg = e12 * sca02 * (_1_0 + 2 * ta02) + e2
                D  = sa0m * (_1_0 - e2 * (_1_0 + sa01 * 2 * sa0)) / (e12 * sa01)  # dD/dsa0
                dD =   -2 * (_1_0 - e2 * sa02 * (_3_0 + 2 * sa0)) / (e12 * sa01**2)
                sa02_ = _1_0 - e2 * sa02
                sa0m_ = sa0m / (_1_0 - e2 * sa0)
                BA = sa0m_ * (self._atanhx1(e2 * sa0m_**2) * e12 - e2 * sa0m) \
                   - sa0m**2 * e2 * (2 + (_1_0 + e2) * sa0) / (e12 * sa02_)  # == B + A
                dAB = 2 * e2 * (2 - e2 * (_1_0 + sa02)) / (e12 * sa02_**2 * sca02)
                u_du = fsum_(s1_qZ *  g,  -D,  g * BA) \
                     / fsum_(s1_qZ * dg, -dD, dg * BA, g * dAB)  # == u/du
                ta0, d = Ta0.fsum2_(-u_du * (sca0 * sca02))
                if abs(d) < tol:
                    break
            else:
                raise AlbersError(iteration=_NUMIT0, txt=_no_(Fmt.convergence(tol)))

        self._txi0  = txi0 =  self._txif(ta0)
        self._scxi0 =         hypot1(txi0)
        self._sxi0  = sxi0 =  txi0 / self._scxi0
        self._m02   = m02  = _1_0 / (_1_0 + (b_a * ta0)**2)
        self._n0    = n0   =  ta0 / hypot1(ta0)
        if polar:
            self._polar = True
            self._nrho0 = self._m0 = _0_0
        else:
            self._m0    = sqrt(m02)  # == nrho0 / E.a
            self._nrho0 = E.a * self._m0  # == E.a * sqrt(m02)
        self._k0_(_1_0 if par1 else (k * sqrt(C / (m02 + n0 * qZ * sxi0))))
        self._lat0 = _Lat(lat0=self._sign * atand(ta0))

    def _s1_qZ_C2(self, ca1, sa1, ta1, ca2, sa2, ta2):
        '''(INTERNAL) Compute C{sm1 / (s / qZ)} and C{C} for .__init__.
        '''
        E = self.datum.ellipsoid
        b_a = E.b_a
        e2  = E.e2

        tb1   =  b_a * ta1
        tb2   =  b_a * ta2
        dtb12 =  b_a * (tb1 + tb2)
        scb12 = _1_0 + tb1**2
        scb22 = _1_0 + tb2**2

        esa1_2 = (_1_0 - e2 * sa1**2) \
               * (_1_0 - e2 * sa2**2)
        esa12  =  _1_0 + e2 * sa1 * sa2

        dsn = _Dsn(ta2, ta1, sa2, sa1)
        axi, bxi, sxi = self._a_b_sxi3((ca1, sa1, ta1, scb12),
                                       (ca2, sa2, ta2, scb22))

        dsxi = (esa12 / esa1_2 + self._Datanhee(sa2, sa1)) * dsn / (_2_0 * self._qx)
        C = fsum_(sxi * dtb12 / dsxi, scb22, scb12) / (_2_0 * scb22 * scb12)

        sa12 = fsum_(sa1, sa2, sa1 * sa2)
        axi *= (_1_0 + e2 * sa12) / (_1_0 + sa12)
        bxi *= e2 * fsum_(sa1, sa2, esa12) / esa1_2 + E.e12 * self._D2atanhee(sa1, sa2)
        s1_qZ = dsn * (axi * self._qZ - bxi) / (_2_0 * dtb12)
        return s1_qZ, C

    def _a_b_sxi3(self, *ca_sa_ta_scb_4s):
        '''(INTERNAL) Sum of C{sm1} terms and C{sin(xi)}s for ._s1_qZ_C2.
        '''
        a = b = s = _0_0
        for ca, sa, ta, scb in ca_sa_ta_scb_4s:
            cxi, sxi, _ = self._cstxif3(ta)
            if sa > 0:
                sa += _1_0
                a  += (cxi / ca)**2 * sa / (_1_0 + sxi)
                b  +=  scb * ca**2 / sa
            else:
                sa =  _1_0 - sa
                a += (_1_0 - sxi) / sa
                b +=  scb * sa
            s += sxi
        return a, b, s

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def equatoradius(self):
        '''Get the geodesic's equatorial radius, semi-axis (C{meter}).
        '''
        return self.datum.ellipsoid.a

    @property_RO
    def iteration(self):
        '''Get the iteration number (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @Property_RO
    def flattening(self):
        '''Get the geodesic's flattening (C{float}).
        '''
        return self.datum.ellipsoid.f

    def forward(self, lat, lon, lon0=0, name=NN):
        '''Convert a geodetic location to east- and northing.

           @arg lat: Latitude of the location (C{degrees}).
           @arg lon: Longitude of the location (C{degrees}).
           @kwarg lon0: Optional central meridian longitude (C{degrees}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Albers7Tuple}C{(x, y, lat, lon, gamma, scale, datum)}.

           @note: The origin latitude is returned by C{property lat0}.  No
                  false easting or northing is added.  The value of B{C{lat}}
                  should be in the range C{[-90..90] degrees}.  The returned
                  values C{x} and C{y} will be large but finite for points
                  projecting to infinity, i.e. one or both of the poles.
        '''
        E = self.datum.ellipsoid
        s = self._sign

        k0    = self._k0
        n0    = self._n0
        nrho0 = self._nrho0
        txi0  = self._txi0

        sa, ca = sincos2d(_Lat(lat=lat) * s)
        ca = max(_EPS__2, ca)
        ta = sa / ca

        _, sxi, txi = self._cstxif3(ta)
        dq = self._qZ * _Dsn(txi, txi0, sxi, self._sxi0) * (txi - txi0)
        drho = -E.a * dq / (sqrt(self._m02 - n0 * dq) + self._m0)

        lon, _ = _diff182(lon0, lon)
        b = radians(lon)

        th = self._k02n0 * b
        sth, cth = sincos2(th)  # XXX sin, cos
        if n0:
            x = sth / n0
            y = ((_1_0 - cth) if cth < 0 else (sth**2 / (_1_0 + cth))) * nrho0 / n0
        else:
            x = self._k02 * b
            y = _0_0
        t = nrho0 + n0 * drho
        x = t * x / k0
        y = s * (y - drho * cth) / k0

        g = degrees360(s * th)
        if t:
            k0 = self._azik(t, ta)
        return Albers7Tuple(x, y, lat, lon, g, k0, self.datum,
                                  name=name or self.name)

    @Property_RO
    def ispolar(self):
        '''Is this projection polar (C{bool})?
        '''
        return self._polar

    isPolar = ispolar  # synonym

    @Property_RO
    def lat0(self):
        '''Get the latitude of the projection origin (C{degrees}).

           This is the latitude of minimum azimuthal scale and
           equals the B{C{lat}} in the 1-parallel L{AlbersEqualArea}
           and lies between B{C{lat1}} and B{C{lat2}} for the
           2-parallel L{AlbersEqualArea2} and L{AlbersEqualArea4}
           projections.
        '''
        return self._lat0

    @Property_RO
    def lat1(self):
        '''Get the latitude of the first parallel (C{degrees}).
        '''
        return self._lat1

    @Property_RO
    def lat2(self):
        '''Get the latitude of the second parallel (C{degrees}).

           @note: The second and first parallel latitudes are the
                  same instance for 1-parallel C{AlbersEqualArea*}
                  projections.
        '''
        return self._lat2

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{equatoradius}.'''
        return self.equatoradius

    def rescale0(self, lat, k=1):
        '''Set the azimuthal scale for this projection.

           @arg lat: Northern latitude (C{degrees}).
           @arg k: Azimuthal scale at latitude B{C{lat}} (C{scalar}).

           @raise AlbersError: Invalid B{C{lat}} or B{C{k}}.

           @note: This allows a I{latitude of conformality} to be specified.
        '''
        k0 = _Ks(k=k) / self.forward(lat, _0_0).scale
        self._update(k0 != self._k0)
        self._k0_(k0)

    def reverse(self, x, y, lon0=0, name=NN, LatLon=None, **LatLon_kwds):
        '''Convert an east- and northing location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg lon0: Optional central meridian longitude (C{degrees}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon}=None}.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Albers7Tuple}C{(x, y, lat, lon, gamma, scale, datum)}.

           @note: The origin latitude is returned by C{property lat0}.  No
                  false easting or northing is added.  The returned value of
                  C{lon} is in the range C{[-180..180] degrees} and C{lat}
                  is in the range C{[-90..90] degrees}.  If the given
                  B{C{x}} or B{C{y}} point is outside the valid projected
                  space the nearest pole is returned.
        '''
        x = Meter(x=x)
        y = Meter(y=y)

        k0    = self._k0
        n0    = self._n0
        k0n0  = self._k0n0
        nrho0 = self._nrho0
        txi0  = self._txi0

        y_ = self._sign * y
        nx = k0n0 * x
        ny = k0n0 * y_
        y1 = nrho0 - ny

        drho = den = nrho0 + hypot(nx, y1)  # 0 implies origin with polar aspect
        if den:
            drho = fsum_(x * nx, -_2_0 * y_ * nrho0, y_ * ny) * k0 / den
        # dsxia = scxi0 * dsxi
        dsxia = -self._scxi0 * (_2_0 * nrho0 + n0 * drho) * drho / self._qZa2
        t = _1_0 - dsxia * (_2_0 * txi0 + dsxia)
        txi = (txi0 + dsxia) / (sqrt(t) if t > _EPS__4 else _EPS__2)

        ta  = self._tanf(txi)
        lat = atand(ta * self._sign)

        th  = atan2(nx, y1)
        lon = degrees(th / self._k02n0 if n0 else x / (y1 * k0))
        if lon0:
            lon += _norm180(lon0)
        lon = _norm180(lon)

        if LatLon is None:
            g = degrees360(self._sign * th)
            if den:
                k0 = self._azik(nrho0 + n0 * drho, ta)
            r = Albers7Tuple(x, y, lat, lon, g, k0, self.datum,
                                   name=name or self.name)
        else:
            kwds = _xkwds(LatLon_kwds, datum=self.datum)
            r = self._xnamed(LatLon(lat, lon, **kwds), name=name)
        return r

    @Property_RO
    def scale0(self):
        '''Get the central scale for the projection (C{float}).

           This is the azimuthal scale on the latitude of origin
           of the projection, see C{property lat0}.
        '''
        return self._k0

    def _atanhee(self, x):
        '''(INTERNAL) Function M{atanhee(x)}, defined as ...
           atanh(      E.e   * x) /       E.e   if f > 0  # oblate
           atan (sqrt(-E.e2) * x) / sqrt(-E.e2) if f < 0  # prolate
           x                                    if f = 0.
        '''
        E = self.datum.ellipsoid
        if E.isOblate:
            x = atanh(x * E.e) / E.e
        elif E.isProlate:
            x = atan2(abs(x) * E.e, copysign0(_1_0, x)) / E.e
        return x

    def _atanhx1(self, x):
        '''(INTERNAL) Function M{atanh(sqrt(x)) / sqrt(x) - 1}.
        '''
        s = abs(x)
        if s < _0_5:  # for typical ...
            # x < E.e^2 = 2 * E.f use ...
            # x / 3 + x^2 / 5 + x^3 / 7 + ...
            y, k = x, 3
            S = Fsum(y / k)
            for _ in range(_TERMS):
                y *= x  # x**n
                k += 2  # 2*n + 1
                s, d = S.fsum2_(y / k)
                if not d:
                    break
        else:
            s =  sqrt(s)
            s = (atanh(s) if x > 0 else atan(s)) / s - _1_0
        return s

    def _azik(self, t, ta):
        '''(INTERNAL) Compute the azimuthal scale C{_Ks(k0=k0)}.
        '''
        E = self.datum.ellipsoid
        return _Ks(k=self._k0 * t * hypot1(E.b_a * ta) / E.a)

    def _cstxif3(self, ta):
        '''(INTERNAL) Get 3-tuple C{(cos, sin, tan)} of M{xi(ta)}.
        '''
        t = self._txif(ta)
        c = _1_0 / hypot1(t)
        s = c * t
        return c, s, t

    def _Datanhee(self, x, y):
        '''(INTERNAL) Function M{Datanhee(x, y)}, defined as
           M{atanhee((x - y) / (1 - E.e^2 * x * y)) / (x - y)}.
        '''
        e2 = self.datum.ellipsoid.e2
        d  = _1_0 - e2 * x * y
        if d:
            d = _1_0 / d
            t =  x - y
            if t:
                d = self._atanhee(t * d) / t
        else:
            raise AlbersError(x=x, y=y, txt=_AlbersBase._Datanhee.__name__)
        return d

    def _D2atanhee(self, x, y):
        '''(INTERNAL) Function M{D2atanhee(x, y)}, defined as
           M{(Datanhee(1, y) - Datanhee(1, x)) / (y - x)}.
        '''
        e2 = self.datum.ellipsoid.e2

        if (abs(x) + abs(y)) * e2 < _0_5:
            e = z = _1_0
            k = 1
            C = Fsum()
            S = Fsum()
            T = Fsum()  # Taylor expansion
            for _ in range(_TERMS):
                T *= y; p = T.fsum_(z); z *= x  # PYCHOK ;
                T *= y; t = T.fsum_(z); z *= x  # PYCHOK ;
                e *= e2
                k += 2
                s, d = S.fsum2_(e * C.fsum_(p, t) / k)
                if not d:
                    break
        elif (_1_0 - x):
            s = (self._Datanhee(_1_0, y) - self._Datanhee(x, y)) / (_1_0 - x)
        else:
            raise AlbersError(x=x, y=y, txt=_AlbersBase._D2atanhee.__name__)
        return s

    def _k0_(self, k):
        '''(INTERNAL) Set C{._k0}, C{._k02}, etc.
        '''
        self._k0    = k  = _Ks(k0=k)
        self._k02   = k2 =  k**2
        self._k0n0  = k  * self._n0
        self._k02n0 = k2 * self._n0

    def _tanf(self, txi):  # called by .Ellipsoid.auxAuthalic
        '''(INTERNAL) Function M{tan-phi from tan-xi}.
        '''
        tol = _tol(_TOL, txi)

        e2 = self.datum.ellipsoid.e2
        qx = self._qx

        ta = txi
        Ta = Fsum(ta)
        for self._iteration in range(1, _NUMIT):  # max 2, mean 1.99
            # dtxi/dta = (scxi / sca)^3 * 2 * (1 - e^2) / (qZ * (1 - e^2 * sa^2)^2)
            ta2  = ta**2
            sca2 = ta2 + _1_0
            txia = self._txif(ta)
            s3qx = sqrt3(sca2 / (_1_0 + txia**2)) * qx
            ta, d = Ta.fsum2_((txi - txia) * s3qx * (_1_0 - e2 * ta2 / sca2)**2)
            if abs(d) < tol:
                return ta
        raise AlbersError(iteration=_NUMIT, txt=_no_(Fmt.convergence(tol)))

    def _txif(self, ta):  # called from .Ellipsoid.auxAuthalic
        '''(INTERNAL) Function M{tan-xi from tan-phi}.
        '''
        E = self.datum.ellipsoid

        ca2 = _1_0 / (_1_0 + ta**2)
        sa  = sqrt(ca2) * abs(ta)  # enforce odd parity

        es1    = E.e2 * sa
        es2m1  = _1_0 - sa * es1
        sp1    = _1_0 + sa
        es1p1  = sp1 / (_1_0 + es1)
        es1m1  = sp1 * (_1_0 - es1)
        es2m1a = es2m1 * E.e12  # e2m
        s = sqrt((ca2 / (es1p1 * es2m1a) + self._atanhee(ca2 / es1m1))
                      * (es1m1 / es2m1a  + self._atanhee(es1p1)))
        t = (sa / es2m1 + self._atanhee(sa)) / s
        if ta < 0:
            t = neg(t)
        return t


class AlbersEqualArea(_AlbersBase):
    '''An Albers equal-area (authalic) projection with a single standard parallel.

       @see: L{AlbersEqualArea2} and L{AlbersEqualArea4}.
    '''
    def __init__(self, lat, k0=1, datum=_WGS84, name=NN):
        '''New L{AlbersEqualArea} projection.

           @arg lat: Standard parallel (C{degrees}).
           @kwarg k0: Azimuthal scale on the standard parallel (C{scalar}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AlbertError: Invalid B{C{lat}}, B{C{k0}} or no convergence.
        '''
        self._lat1 = self._lat2 = lat = _Lat(lat1=lat)
        args = tuple(sincos2d(lat)) * 2 + (_Ks(k0=k0), datum, name)
        _AlbersBase.__init__(self, *args)


class AlbersEqualArea2(_AlbersBase):
    '''An Albers equal-area (authalic) projection with two standard parallels.

       @see: L{AlbersEqualArea} and L{AlbersEqualArea4}.
    '''
    def __init__(self, lat1, lat2, k1=1, datum=_WGS84, name=NN):
        '''New L{AlbersEqualArea2} projection.

           @arg lat1: First standard parallel (C{degrees}).
           @arg lat2: Second standard parallel (C{degrees}).
           @kwarg k1: Azimuthal scale on the standard parallels (C{scalar}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AlbertError: Invalid B{C{lat1}}m B{C{lat2}}, B{C{k1}}
                               or no convergence.
        '''
        self._lat1, self._lat2 = lats = _Lat(lat1=lat1), _Lat(lat2=lat2)
        args = tuple(sincos2d(*lats)) + (_Ks(k1=k1), datum, name)
        _AlbersBase.__init__(self, *args)


class AlbersEqualArea4(_AlbersBase):
    '''An Albers equal-area (authalic) projection specified by the C{sin}
       and C{cos} of both standard parallels.

       @see: L{AlbersEqualArea} and L{AlbersEqualArea2}.
    '''
    def __init__(self, slat1, clat1, slat2, clat2, k1=1, datum=_WGS84, name=NN):
        '''New L{AlbersEqualArea4} projection.

           @arg slat1: Sine of first standard parallel (C{scalar}).
           @arg clat1: Cosine of first standard parallel (non-negative C{scalar}).
           @arg slat2: Sine of second standard parallel (C{scalar}).
           @arg clat2: Cosine of second standard parallel (non-negative C{scalar}).
           @kwarg k1: Azimuthal scale on the standard parallels (C{scalar}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AlbertError: Negative B{C{clat1}} or B{C{clat2}}, B{C{slat1}}
                               and B{C{slat2}} have opposite signs (hemispheres),
                               invalid B{C{k1}} or no convergence.
        '''
        def _Lat_s_c3(name, s, c):
            L =  Lat_(atan2d(s, c), name=name, Error=AlbersError)
            r = _1_0 / Float_(hypot(s, c), name=name, Error=AlbersError)
            return L, s * r, c * r

        self._lat1, sa1, ca1 = _Lat_s_c3(_lat1_, slat1, clat1)
        self._lat2, sa2, ca2 = _Lat_s_c3(_lat2_, slat2, clat2)
        _AlbersBase.__init__(self, sa1, ca1, sa2, ca2, _Ks(k1=k1), datum, name)


class AlbersEqualAreaCylindrical(_AlbersBase):
    '''An L{AlbersEqualArea} projection at C{lat=0} and C{k0=1} degenerating
       to the cylindrical-equal-area projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=_0_0)

    def __init__(self, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaCylindrical} projection.

           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _AlbersBase.__init__(self, _0_0, _1_0, _0_0, _1_0, _1_0, datum, name)


class AlbersEqualAreaNorth(_AlbersBase):
    '''An azimuthal L{AlbersEqualArea} projection at C{lat=90} and C{k0=1}
       degenerating to the L{azimuthal} L{LambertEqualArea} projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=_90_0)

    def __init__(self, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaNorth} projection.

           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _AlbersBase.__init__(self, _1_0, _0_0, _1_0, _0_0, _1_0, datum, name)


class AlbersEqualAreaSouth(_AlbersBase):
    '''An azimuthal L{AlbersEqualArea} projection at C{lat=-90} and C{k0=1}
       degenerating to the L{azimuthal} L{LambertEqualArea} projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=-_90_0)

    def __init__(self, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaSouth} projection.

           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _AlbersBase.__init__(self, -_1_0, _0_0, -_1_0, _0_0, _1_0, datum, name)


class Albers7Tuple(_NamedTuple):
    '''7-Tuple C{(x, y, lat, lon, gamma, scale, datum)}, in C{meter},
       C{meter}, C{degrees90}, C{degrees180}, C{degrees360}, C{scalar} and
       C{Datum} where C{(x, y)} is the projected, C{(lat, lon)} the geodetic
       location, C{gamma} the meridian convergence at point, the bearing of
       the y-axis measured clockwise from true North and C{scale} is the
       azimuthal scale of the projection at point.  The radial scale is
       the reciprocal C{1 / scale}.
    '''
    _Names_ = (_x_,   _y_,   _lat_, _lon_, _gamma_,  _scale_, _datum_)
    _Units_ = ( Meter, Meter, Lat,   Lon,   Bearing, _Pass,   _Pass)


__all__ += _ALL_DOCS(_AlbersBase)

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
