
# -*- coding: utf-8 -*-

u'''Albers Equal-Area projections.

Classes L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth}, L{AlbersEqualAreaSouth}
and L{AlbersError}, transcoded from I{Charles Karney}'s C++ class U{AlbersEqualArea
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AlbersEqualArea.html>}.

See also I{Albers Equal-Area Conic Projection} in U{John P. Snyder, "Map Projections
-- A Working Manual", 1987<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 98-106
and the Albers Conical Equal-Area examples on pp 291-294.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import neg, neg_
from pygeodesy.constants import EPS0, EPS02, _EPSqrt as _TOL, \
                               _0_0, _0_5, _1_0, _N_1_0, _2_0, \
                               _N_2_0, _4_0, _6_0, _90_0, _N_90_0
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _ValueError, _xkwds
from pygeodesy.fmath import hypot, hypot1, sqrt3
from pygeodesy.fsums import Fsum, fsum1f_
from pygeodesy.interns import NN, _COMMASPACE_, _datum_, _gamma_, _k0_, \
                             _lat_, _lat1_, _lat2_, _lon_, _name_, _not_, \
                             _negative_, _scale_, _SPACE_, _x_, _y_
from pygeodesy.karney import _diff182, _norm180, _signBit
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass
from pygeodesy.props import deprecated_Property_RO, Property_RO, _update_all
from pygeodesy.streprs import Fmt, strs, unstr
from pygeodesy.units import Bearing, Float_, Lat, Lat_, Lon, Meter, Scalar_
from pygeodesy.utily import atand, atan2d, degrees360, sincos2, \
                            sincos2d, sincos2d_

from math import atan, atan2, atanh, degrees, fabs, radians, sqrt

__all__ = _ALL_LAZY.albers
__version__ = '23.09.07'

_k1_    = 'k1'
_NUMIT  =   8  # XXX 4?
_NUMIT0 =  41  # XXX 21?
_TERMS  =  21  # XXX 16?
_TOL0   =  sqrt3(_TOL)


def _ct2(s, c):
    '''(INTERNAL) Avoid singularities at poles.
    '''
    c = max(EPS0, c)
    return c, (s / c)


def _Fsum1(a, b, c=_0_0):  # floats=True
    '''(INTERNAL) C{Fsum} 1-primed.
    '''
    S = Fsum()
    S._facc_(_1_0, a, b, c, _N_1_0, up=False)
    return S


def _Ks(**name_k):
    '''(INTERNAL) Scale C{B{k} >= EPS0}.
    '''
    return Scalar_(Error=AlbersError, low=EPS0, **name_k)  # > 0


def _Lat(*lat, **Error_name_lat):
    '''(INTERNAL) Latitude C{-90 <= B{lat} <= 90}.
    '''
    kwds = _xkwds(Error_name_lat, Error=AlbersError)
    return Lat_(*lat, **kwds)


def _qZx(albs):
    '''(INTERNAL) Set C{albs._qZ} and C{albs._qx}.
    '''
    E = albs._datum.ellipsoid  # _AlbersBase
    albs._qZ = qZ =  _1_0 + E.e21 * _atanheE(_1_0, E)
    albs._qx = qZ / (_2_0 * E.e21)
    return qZ


class AlbersError(_ValueError):
    '''An L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
       L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth},
       L{AlbersEqualAreaSouth} or L{Albers7Tuple} issue.
    '''
    pass


class _AlbersBase(_NamedBase):
    '''(INTERNAL) Base class for C{AlbersEqualArea...} projections.

       @see: I{Karney}'s C++ class U{AlbersEqualArea<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1AlbersEqualArea.html>}, method C{Init}.
    '''
    _datum  = _WGS84
    _k      =  NN    # or _k0_ or _k1_
    _k0     = _Ks(k0=_1_0)
#   _k0n0   =  None  # (INTERNAL) k0 * no
    _k02    = _1_0   # (INTERNAL) k0**2
#   _k02n0  =  None  # (INTERNAL) k02 * n0
#   _lat0   =  None  # lat origin
    _lat1   =  None  # let 1st parallel
    _lat2   =  None  # lat 2nd parallel
    _m0     = _0_0   # if polar else sqrt(m02)
#   _m02    =  None  # (INTERNAL) cached
#   _n0     =  None  # (INTERNAL) cached
    _nrho0  = _0_0   # if polar else m0 * E.a
    _polar  =  False
    _qx     =  None  # (INTERNAL) see _qZx
    _qZ     =  None  # (INTERNAL) see _qZx
#   _scxi0_ =  None  # (INTERNAL) sec(xi) / (qZ * E.a2)
    _sign   = +1
#   _sxi0   =  None  # (INTERNAL) sin(xi)
#   _txi0   =  None  # (INTERNAL) tan(xi)

    def __init__(self, sa1, ca1, sa2, ca2, k, datum, name):
        '''(INTERNAL) New C{AlbersEqualArea...} instance.
        '''
        qZ = self._qZ
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)
            qZ = _qZx(self)
        elif qZ is None:
            qZ = _qZx(_AlbersBase)
        if name:
            self.name = name

        E = self.ellipsoid
        c = min(ca1, ca2)
        if _signBit(c):
            raise AlbersError(clat1=ca1, clat2=ca2, txt=_negative_)
        polar = c < EPS0  # == 0

        # determine hemisphere of tangent latitude
        if sa1 < 0:  # and sa2 < 0:
            self._sign = -1
            # internally, tangent latitude positive
            sa1, sa2 = neg_(sa1, sa2)
        if sa1 > sa2:  # make phi1 < phi2
            sa1, sa2 = sa2, sa1
            ca1, ca2 = ca2, ca1
        if sa1 < 0:  # or sa2 < 0:
            raise AlbersError(slat1=sa1, slat2=sa2, txt=_negative_)
        ca1, ta1 = _ct2(sa1, ca1)
        ca2, ta2 = _ct2(sa2, ca2)

        par1 = fabs(ta1 - ta2) < EPS02  # ta1 == ta2
        if par1 or polar:
            ta0, C = ta2, _1_0
        else:
            ta0, C = self._ta0C2(ca1, sa1, ta1, ca2, sa2, ta2)

        self._lat0 = _Lat(lat0=self._sign * atand(ta0))
        self._m02  =  m02 = _1_x21(E.f1 * ta0)
        self._n0   =  n0  =  ta0 / hypot1(ta0)
        if polar:
            self._polar = True
#           self._nrho0 = self._m0 = _0_0
        else:  # m0 = nrho0 / E.a
            self._m0    = sqrt(m02)
            self._nrho0 = self._m0 * E.a
        t = self._txi0  = self._txif(ta0)
        h = hypot1(t)
        s = self._sxi0  = t / h
        if par1:
            self._k0n0  = self._k02n0 = n0
        else:
            self._k0s(k * sqrt(C / (m02 + n0 * qZ * s)))
        self._scxi0_ = h / (qZ * E.a2)

    def _a_b_sxi3(self, *ca_sa_ta_scb_4s):
        '''(INTERNAL) Sum of C{sm1} terms and C{sin(xi)}s for ._ta0C2.
        '''
        _1 = _1_0
        a  =  b = s = _0_0
        for ca, sa, ta, scb in ca_sa_ta_scb_4s:
            cxi, sxi, _ = self._cstxif3(ta)
            if sa > 0:
                sa += _1
                a  += (cxi / ca)**2 * sa / (sxi + _1)
                b  +=  scb * ca**2  / sa
            else:
                sa =  _1 - sa
                a += (_1 - sxi) / sa
                b +=  scb * sa
            s += sxi
        return a, b, s

    def _azik(self, t, ta):
        '''(INTERNAL) Compute the azimuthal scale C{_Ks(k=k)}.
        '''
        E = self.ellipsoid
        return _Ks(k=hypot1(E.b_a * ta) * self._k0 * t / E.a)

    def _cstxif3(self, ta):
        '''(INTERNAL) Get 3-tuple C{(cos, sin, tan)} of M{xi(ta)}.
        '''
        t =  self._txif(ta)
        c = _1_0 / hypot1(t)
        s =  c * t
        return c, s, t

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def ellipsoid(self):
        '''Get the datum's ellipsoid (L{Ellipsoid}).
        '''
        return self.datum.ellipsoid

    @Property_RO
    def equatoradius(self):
        '''Get the C{ellipsoid}'s equatorial radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    a = equatoradius

    @Property_RO
    def flattening(self):
        '''Get the C{ellipsoid}'s flattening (C{scalar}).
        '''
        return self.ellipsoid.f

    f = flattening

    def forward(self, lat, lon, lon0=0, name=NN):
        '''Convert a geodetic location to east- and northing.

           @arg lat: Latitude of the location (C{degrees}).
           @arg lon: Longitude of the location (C{degrees}).
           @kwarg lon0: Optional central meridian longitude (C{degrees}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Albers7Tuple}C{(x, y, lat, lon, gamma, scale, datum)},
                    with C{lon} offset by B{C{lon0}} and reduced C{[-180,180]}.

           @note: The origin latitude is returned by C{property lat0}.  No
                  false easting or northing is added.  The value of B{C{lat}}
                  should be in the range C{[-90..90] degrees}.  The returned
                  values C{x} and C{y} will be large but finite for points
                  projecting to infinity, i.e. one or both of the poles.
        '''
        a = self.ellipsoid.a
        s = self._sign

        k0    = self._k0
        n0    = self._n0
        nrho0 = self._nrho0
        txi0  = self._txi0

        _, ta = _ct2(*sincos2d(s * _Lat(lat=lat)))

        _, sxi, txi = self._cstxif3(ta)
        dq = _Dsn(txi,  txi0, sxi, self._sxi0) * \
                 (txi - txi0) * self._qZ
        drho = a * dq / (sqrt(self._m02 - n0 * dq) + self._m0)

        lon, _ = _diff182(lon0, lon)
        x  = radians(lon)
        th = self._k02n0 * x
        sth, cth = sincos2(th)  # XXX sin, cos
        if n0:
            x  = sth / n0
            y  = (_1_0 - cth) if cth < 0 else (sth**2 / (cth + _1_0))
            y *= nrho0 / n0
        else:
            x *= self._k02
            y = _0_0
        t = nrho0 - n0 * drho
        x = t *  x / k0
        y = s * (y + drho * cth) / k0

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

    def _k0s(self, k0):
        '''(INTERNAL) Set C{._k0}, C{._k02}, etc.
        '''
        self._k0    = k0  = _Ks(k0=k0)
        self._k02   = k02 =  k0**2
        self._k0n0  = k0  *  self._n0
        self._k02n0 = k02 *  self._n0

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

    def rescale0(self, lat, k=1):  # PYCHOK no cover
        '''Set the azimuthal scale for this projection.

           @arg lat: Northern latitude (C{degrees}).
           @arg k: Azimuthal scale at latitude B{C{lat}} (C{scalar}).

           @raise AlbersError: Invalid B{C{lat}} or B{C{k}}.

           @note: This allows a I{latitude of conformality} to be specified.
        '''
        k0 = _Ks(k=k) / self.forward(lat, _0_0).scale
        if self._k0 != k0:
            _update_all(self)
            self._k0s(k0)

    def reverse(self, x, y, lon0=0, name=NN, LatLon=None, **LatLon_kwds):
        '''Convert an east- and northing location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg lon0: Optional central meridian longitude (C{degrees}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Albers7Tuple}C{(x, y, lat, lon, gamma, scale, datum)}.

           @note: The origin latitude is returned by C{property lat0}.  No
                  false easting or northing is added.  The returned value of
                  C{lon} is in the range C{[-180..180] degrees} and C{lat}
                  is in the range C{[-90..90] degrees}.  If the given
                  B{C{x}} or B{C{y}} point is outside the valid projected
                  space the nearest pole is returned.
        '''
        k0   = self._k0
        n0   = self._n0
        k0n0 = self._k0n0
        s    = self._sign
        txi  = self._txi0

        x  = Meter(x=x)
        nx = k0n0 * x
        y  = Meter(y=y)
        y_ = s    * y
        ny = k0n0 * y_
        t  = nrho0 = self._nrho0
        y1 = nrho0 - ny

        den = hypot(nx, y1) + nrho0  # 0 implies origin with polar aspect
        if den:
            drho = _Fsum1(x * nx, y_ * nrho0 * _N_2_0, y_ * ny).fover(den / k0)
            # dsxia = scxi0 * dsxi
            t  +=  drho  * n0
            d_  = (nrho0 + t) * drho * self._scxi0_  # / (qZ * E.a2)
            d_2 = (txi * _2_0 - d_) * d_ + _1_0
            txi = (txi - d_) / (sqrt(d_2) if d_2 > EPS02 else EPS0)

        ta  = self._tanf(txi)
        lat = atand(s * ta)

        th  = atan2(nx, y1)
        lon = degrees((th / self._k02n0) if n0 else (x / (y1 * k0)))
        if lon0:
            lon += _norm180(lon0)
        lon = _norm180(lon)

        n = name or self.name
        if LatLon is None:
            g = degrees360(s * th)
            if den:
                k0 = self._azik(t, ta)
            r = Albers7Tuple(x, y, lat, lon, g, k0, self.datum, name=n)
        else:  # PYCHOK no cover
            kwds = _xkwds(LatLon_kwds, datum=self.datum, name=n)
            r = LatLon(lat, lon, **kwds)
        return r

    @Property_RO
    def scale0(self):
        '''Get the central scale for the projection (C{float}).

           This is the azimuthal scale on the latitude of origin
           of the projection, see C{property lat0}.
        '''
        return self._k0

    def _ta0(self, s1_qZ, ta0, E):
        '''(INTERNAL) Refine C{ta0} for C{._ta0C2}.
        '''
        e2     =  E.e2
        e21    =  E.e21
        e22    =  E.e22  # == e2 / e21
        tol    = _tol(_TOL0, ta0)
        _Ta02  =  Fsum(ta0).fsum2_
        _fabs  =  fabs
        _fsum1 =  fsum1f_
        _sqrt  =  sqrt
        _1, _2 = _1_0, _2_0
        _4, _6 = _4_0, _6_0
        for self._iteration in range(1, _NUMIT0):  # 4 trips
            ta02  =  ta0**2
            sca02 =  ta02 + _1
            sca0  = _sqrt(sca02)
            sa0   =  ta0 / sca0
            sa01  =  sa0 + _1
            sa02  =  sa0**2
            # sa0m = 1 - sa0 = 1 / (sec(a0) * (tan(a0) + sec(a0)))
            sa0m  = _1 / (sca0 * (ta0 + sca0))  # scb0^2 * sa0
            sa0m1 =  sa0m / (_1 - e2 * sa0)
            sa021 =          _1 - e2 * sa02

            g  = (_1 + ta02 * e21) * sa0
            dg = (_1 + ta02 * _2)  * sca02 * e21 + e2
            D  = (_1 - (_1 + sa0 * _2  * sa01) * e2) * sa0m / (e21 * sa01)  # dD/dsa0
            dD = (_2 - (_6 + sa0 * _4) * sa02  * e2)       /  (e21 * sa01**2)
            BA = (_atanh1(e2 * sa0m1**2) * e21 - e2 * sa0m) * sa0m1 \
               - (_2 + (_1 + e2) * sa0) * sa0m**2 * e22 / sa021  # B + A
            d  = (_4 - (_1 + sa02) * e2 * _2) * e22    / (sa021**2 * sca02)  # dAB
            u  = _fsum1(s1_qZ *  g, -D,  g * BA)
            du = _fsum1(s1_qZ * dg, dD, dg * BA, g * d)
            ta0, d = _Ta02(-u / du * (sca0 * sca02))
            if _fabs(d) < tol:
                return ta0
        raise AlbersError(Fmt.no_convergence(d, tol), txt=repr(self))

    def _ta0C2(self, ca1, sa1, ta1, ca2, sa2, ta2):
        '''(INTERNAL) Compute C{ta0} and C{C} for C{.__init__}.
        '''
        E      =  self.ellipsoid
        f1, e2 =  E.f1, E.e2
        _1     = _1_0

        tb1    = f1 *  ta1
        tb2    = f1 *  ta2
        dtb12  = f1 * (tb1 + tb2)
        scb12  = _1 +  tb1**2
        scb22  = _1 +  tb2**2

        dsn_2  = _Dsn(ta2, ta1, sa2, sa1) * _0_5
        sa12   =  sa1 * sa2

        esa1_2 = (_1 - e2 * sa1**2) \
               * (_1 - e2 * sa2**2)
        esa12  =  _1 + e2 * sa12

        axi, bxi, sxi = self._a_b_sxi3((ca1, sa1, ta1, scb12),
                                       (ca2, sa2, ta2, scb22))

        dsxi = ((esa12 / esa1_2) + _DatanheE(sa2, sa1, E)) * dsn_2 / self._qx
        C = _Fsum1(sxi * dtb12 / dsxi, scb22, scb12).fover(scb22 * scb12 * _2_0)

        sa12  =  fsum1f_(sa1, sa2, sa12)
        axi  *= (sa12 * e2 + _1) / (sa12 + _1)
        bxi  *= _Fsum1(sa1, sa2, esa12).fover(esa1_2) * e2 + _D2atanheE(sa1, sa2, E) * E.e21
        s1_qZ = (axi * self._qZ - bxi) * dsn_2 / dtb12
        ta0   =  self._ta0(s1_qZ, (ta1 + ta2) * _0_5, E)
        return ta0, C

    def _tanf(self, txi):  # in .Ellipsoid.auxAuthalic
        '''(INTERNAL) Function M{tan-phi from tan-xi}.
        '''
        tol = _tol(_TOL, txi)
        e2  =  self.ellipsoid.e2
        qx  =  self._qx

        ta     =  txi
        _Ta2   =  Fsum(ta).fsum2_
        _fabs  =  fabs
        _sqrt3 =  sqrt3
        _txif  =  self._txif
        _1     = _1_0
        for self._iteration in range(1, _NUMIT):  # max 2, mean 1.99
            # dtxi / dta = (scxi / sca)^3 * 2 * (1 - e^2)
            #            / (qZ * (1 - e^2 * sa^2)^2)
            ta2  =  ta**2
            sca2 = _1 + ta2
            txia = _txif(ta)
            s3qx = _sqrt3(sca2 / (txia**2 + _1)) * qx  # * _1_x21(txia)
            eta2 = (_1 - e2 * ta2 / sca2)**2
            ta, d = _Ta2((txi - txia) * s3qx * eta2)
            if _fabs(d) < tol:
                return ta
        raise AlbersError(Fmt.no_convergence(d, tol), txt=repr(self))

    def toRepr(self, prec=6, **unused):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).

           @return: This projection as C{"<classname>(lat1, lat2, ...)"}
                    (C{str}).
        '''
        t = self.toStr(prec=prec, sep=_COMMASPACE_)
        return Fmt.PAREN(self.classname, t)

    def toStr(self, prec=6, sep=_SPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg sep: Separator to join (C{str}).

           @return: This projection as C{"lat1 lat2"} (C{str}).
        '''
        k =  self._k
        t = (self.lat1, self.lat2, self._k0) if k is _k1_ else (
            (self.lat1,            self._k0) if k is _k0_ else
            (self.lat1,))
        t = strs(t, prec=prec)
        if k:
            t  =  t[:-1] + (Fmt.EQUAL(k, t[-1]),)
        if self.datum != _WGS84:
            t += (Fmt.EQUAL(_datum_, self.datum),)
        if self.name:
            t += (Fmt.EQUAL(_name_, repr(self.name)),)
        return t if sep is None else sep.join(t)

    def _txif(self, ta):  # in .Ellipsoid.auxAuthalic
        '''(INTERNAL) Function M{tan-xi from tan-phi}.
        '''
        E   =  self.ellipsoid
        _1  = _1_0

        ca2 = _1_x21(ta)
        sa  =  sqrt(ca2) * fabs(ta)  # enforce odd parity
        sa1 = _1 + sa

        es1    =  sa  * E.e2
        es1m1  =  sa1 * (_1 - es1)
        es1p1  =  sa1 / (_1 + es1)
        es2m1  = _1   -  sa * es1
        es2m1a =  es2m1 * E.e21  # e2m
        s =  sqrt((ca2 / (es1p1 * es2m1a) + _atanheE(ca2 / es1m1, E))
                       * (es1m1 / es2m1a  + _atanheE(es1p1, E)))
        t = _Fsum1(sa / es2m1, _atanheE(sa, E)).fover(s)
        return neg(t) if ta < 0 else t


class AlbersEqualArea(_AlbersBase):
    '''An Albers equal-area (authalic) projection with a single standard parallel.

       @see: L{AlbersEqualArea2} and L{AlbersEqualArea4}.
    '''
    _k = _k0_

    def __init__(self, lat, k0=1, datum=_WGS84, name=NN):
        '''New L{AlbersEqualArea} projection.

           @arg lat: Standard parallel (C{degrees}).
           @kwarg k0: Azimuthal scale on the standard parallel (C{scalar}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AlbersError: Invalid B{C{lat}}, B{C{k0}} or no convergence.
        '''
        self._lat1 = self._lat2 = lat = _Lat(lat1=lat)
        args = tuple(sincos2d(lat)) * 2 + (_Ks(k0=k0), datum, name)
        _AlbersBase.__init__(self, *args)


class AlbersEqualArea2(_AlbersBase):
    '''An Albers equal-area (authalic) projection with two standard parallels.

       @see: L{AlbersEqualArea} and L{AlbersEqualArea4}.
    '''
    _k = _k1_

    def __init__(self, lat1, lat2, k1=1, datum=_WGS84, name=NN):
        '''New L{AlbersEqualArea2} projection.

           @arg lat1: First standard parallel (C{degrees}).
           @arg lat2: Second standard parallel (C{degrees}).
           @kwarg k1: Azimuthal scale on the standard parallels (C{scalar}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AlbersError: Invalid B{C{lat1}}m B{C{lat2}}, B{C{k1}}
                               or no convergence.
        '''
        self._lat1, self._lat2 = lats = _Lat(lat1=lat1), _Lat(lat2=lat2)
        args = tuple(sincos2d_(*lats)) + (_Ks(k1=k1), datum, name)
        _AlbersBase.__init__(self, *args)


class AlbersEqualArea4(_AlbersBase):
    '''An Albers equal-area (authalic) projection specified by the C{sin}
       and C{cos} of both standard parallels.

       @see: L{AlbersEqualArea} and L{AlbersEqualArea2}.
    '''
    _k = _k1_

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

           @raise AlbersError: Negative B{C{clat1}} or B{C{clat2}}, B{C{slat1}}
                               and B{C{slat2}} have opposite signs (hemispheres),
                               invalid B{C{k1}} or no convergence.
        '''
        def _Lat_s_c3(name, s, c):
            r =  Float_(hypot(s, c), name=name, Error=AlbersError)
            L = _Lat(atan2d(s, c), name=name)
            return L, (s / r), (c / r)

        self._lat1, sa1, ca1 = _Lat_s_c3(_lat1_, slat1, clat1)
        self._lat2, sa2, ca2 = _Lat_s_c3(_lat2_, slat2, clat2)
        _AlbersBase.__init__(self, sa1, ca1, sa2, ca2, _Ks(k1=k1), datum, name)


class AlbersEqualAreaCylindrical(_AlbersBase):
    '''An L{AlbersEqualArea} projection at C{lat=0} and C{k0=1} degenerating
       to the cylindrical-equal-area projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=_0_0)

    def __init__(self, lat=_0_0, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaCylindrical} projection.

           @kwarg lat: Standard parallel (C{0 degrees} I{fixed}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _xlat(lat, _0_0, AlbersEqualAreaCylindrical)
        _AlbersBase.__init__(self, _0_0, _1_0, _0_0, _1_0, 1, datum, name)


class AlbersEqualAreaNorth(_AlbersBase):
    '''An azimuthal L{AlbersEqualArea} projection at C{lat=90} and C{k0=1}
       degenerating to the L{azimuthal} L{LambertEqualArea} projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=_90_0)

    def __init__(self, lat=_90_0, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaNorth} projection.

           @kwarg lat: Standard parallel (C{90 degrees} I{fixed}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _xlat(lat, _90_0, AlbersEqualAreaNorth)
        _AlbersBase.__init__(self, _1_0, _0_0, _1_0, _0_0, 1, datum, name)


class AlbersEqualAreaSouth(_AlbersBase):
    '''An azimuthal L{AlbersEqualArea} projection at C{lat=-90} and C{k0=1}
       degenerating to the L{azimuthal} L{LambertEqualArea} projection.
    '''
    _lat1 = _lat2 = _Lat(lat1=_N_90_0)

    def __init__(self, lat=_N_90_0, datum=_WGS84, name=NN):
        '''New L{AlbersEqualAreaSouth} projection.

           @kwarg lat: Standard parallel (C{-90 degrees} I{fixed}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name for the projection (C{str}).
        '''
        _xlat(lat, _N_90_0, AlbersEqualAreaSouth)
        _AlbersBase.__init__(self, _N_1_0, _0_0, _N_1_0, _0_0, 1, datum, name)


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


def _atanh1(x):
    '''(INTERNAL) Function M{atanh(sqrt(x)) / sqrt(x) - 1}.
    '''
    s = fabs(x)
    if 0 < s < _0_5:  # for typical ...
        # x < E.e^2 == 2 * E.f use ...
        # x / 3 + x^2 / 5 + x^3 / 7 + ...
        y, k = x, 3
        _S2  = Fsum(y / k).fsum2_
        for _ in range(_TERMS):  # 9 terms
            y *= x  # x**n
            k += 2  # 2*n + 1
            s, d = _S2(y / k)
            if not d:
                break
    elif s:
        s =  sqrt(s)
        s = (atanh(s) if x > 0 else atan(s)) / s - _1_0
    return s


def _atanheE(x, E):  # see Ellipsoid._es_atanh, .AuxLat._atanhee
    '''(INTERNAL) Function M{atanhee(x)}, defined as ...
       atanh(      E.e   * x) /       E.e   if f > 0  # oblate
       atan (sqrt(-E.e2) * x) / sqrt(-E.e2) if f < 0  # prolate
       x                                    if f = 0.
    '''
    e = E.e
    if E.f > 0:  # .isOblate
        x =  atanh( x * e) / e
    elif E.f < 0:  # .isProlate
        x = (atan2(-x * e, _N_1_0) if x < 0 else
             atan2( x * e,   _1_0)) / e
    return x


def _DatanheE(x, y, E):  # see .rhumbx._DeatanhE
    '''(INTERNAL) Function M{Datanhee(x, y)}, defined as
       M{atanhee((x - y) / (1 - E.e^2 * x * y)) / (x - y)}.
    '''
    e = _1_0 - E.e2 * x * y
    if e:
        d = x - y
        e = (_atanheE(d / e, E) / d) if d else (_1_0 / e)
    return e


def _D2atanheE(x, y, E):
    '''(INTERNAL) Function M{D2atanhee(x, y)}, defined as
       M{(Datanhee(1, y) - Datanhee(1, x)) / (y - x)}.
    '''
    s, e2 = _0_0, E.e2
    if e2:
        if ((fabs(x) + fabs(y)) * e2) < _0_5:
            e = z = _1_0
            k = 1
            T   = Fsum()  # Taylor expansion
            _T  = T.fsum_
            _C  = Fsum().fsum_
            _S2 = Fsum().fsum2_
            for _ in range(_TERMS):  # 15 terms
                T *= y; p = _T(z); z *= x  # PYCHOK ;
                T *= y; t = _T(z); z *= x  # PYCHOK ;
                e *= e2
                k += 2
                s, d = _S2(_C(p, t) * e / k)
                if not d:
                    break
        else:  # PYCHOK no cover
            s = _1_0 - x
            if s:
                s = (_DatanheE(_1_0, y, E) - _DatanheE(x, y, E)) / s
    return s


def _Dsn(x, y, sx, sy):
    '''(INTERNAL) Divided differences, defined as M{Df(x, y) = (f(x) - f(y)) / (x - y)}
       with M{sn(x) = x / sqrt(1 + x^2)}: M{Dsn(x, y) = (x + y) / ((sn(x) + sn(y)) *
       (1 + x^2) * (1 + y^2))}.

       @see: U{W. M. Kahan and R. J. Fateman, "Sympbolic Computation of Divided
             Differences"<https://People.EECS.Berkeley.EDU/~fateman/papers/divdiff.pdf>},
             U{ACM SIGSAM Bulletin 33(2), 7-28 (1999)<https://DOI.org/10.1145/334714.334716>}
             and U{AlbersEqualArea.hpp
             <https://GeographicLib.SourceForge.io/C++/doc/AlbersEqualArea_8hpp_source.html>}.
    '''
    # sx = x / hypot1(x)
    d, t = _1_0, (x * y)
    if t > 0:
        s = sx + sy
        if s:
            t = sx * sy / t
            d = t**2 * (x + y) / s
    elif x != y:
        d = (sx - sy) / (x - y)
    return d


def _tol(tol, x):
    '''(INTERNAL) Converge tolerance.
    '''
    return tol * max(_1_0, fabs(x))


def _1_x21(x):
    '''(INTERNAL) Return M{1 / (x**2 + 1)}.
    '''
    return _1_0 / (x**2 + _1_0)


def _xlat(lat, f, where):
    '''(INTERNAL) check fixed C{lat}.
    '''
    if lat is not f and _Lat(lat=lat) != f:
        t = unstr(where.__name__, lat=lat)
        raise AlbersError(t, txt=_not_(f))


__all__ += _ALL_DOCS(_AlbersBase)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
