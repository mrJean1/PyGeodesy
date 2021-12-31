
# -*- coding: utf-8 -*-

u'''Lambert Conformal Conic (LCC) projection.

Lambert conformal conic projection for 1- or 2-Standard Parallels classes L{Conic}, L{Conics} registry, L{LCCError}
and position class L{Lcc}.

See U{LCC<https://WikiPedia.org/wiki/Lambert_conformal_conic_projection>}, U{Lambert
Conformal Conic to Geographic Transformation Formulae
<https://www.Linz.govt.NZ/data/geodetic-system/coordinate-conversion/projection-conversions/lambert-conformal-conic-geographic>},
U{Lambert Conformal Conic Projection<https://MathWorld.Wolfram.com/LambertConformalConicProjection.html>}
and John P. Snyder U{'Map Projections - A Working Manual'<https://Pubs.USGS.gov/pp/1395/report.pdf>}, 1987, pp 107-109.

@var Conics.Be08Lb: Conic(name='Be08Lb', lat0=50.797815, lon0=4.35921583, par1=49.8333339, par2=51.1666672, E0=649328, N0=665262, k0=1, SP=2, datum=Datum(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84),
@var Conics.Be72Lb: Conic(name='Be72Lb', lat0=90, lon0=4.3674867, par1=49.8333339, par2=51.1666672, E0=150000.013, N0=5400088.438, k0=1, SP=2, datum=Datum(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83),
@var Conics.Fr93Lb: Conic(name='Fr93Lb', lat0=46.5, lon0=3, par1=49, par2=44, E0=700000, N0=6600000, k0=1, SP=2, datum=Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84),
@var Conics.MaNLb: Conic(name='MaNLb', lat0=33.3, lon0=-5.4, par1=31.73, par2=34.87, E0=500000, N0=300000, k0=1, SP=2, datum=Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF),
@var Conics.MxLb: Conic(name='MxLb', lat0=12, lon0=-102, par1=17.5, par2=29.5, E0=2500000, N0=0, k0=1, SP=2, datum=Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84),
@var Conics.PyT_Lb: Conic(name='PyT_Lb', lat0=46.8, lon0=2.33722917, par1=45.8989389, par2=47.6960144, E0=600000, N0=200000, k0=1, SP=2, datum=Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF),
@var Conics.USA_Lb: Conic(name='USA_Lb', lat0=23, lon0=-96, par1=33, par2=45, E0=0, N0=0, k0=1, SP=2, datum=Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84),
@var Conics.WRF_Lb: Conic(name='WRF_Lb', lat0=40, lon0=-97, par1=33, par2=45, E0=0, N0=0, k0=1, SP=2, datum=Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, _xinstanceof, _xsubclassof
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.datums import Datums, _ellipsoidal_datum
from pygeodesy.errors import _IsnotError, _ValueError
from pygeodesy.fmath import hypot
from pygeodesy.interns import EPS, EPS02, NN, PI_2, _COMMASPACE_, _ellipsoidal_, \
                             _float as _F, _GRS80_, _k0_, _lat0_, _lon0_, _m_, \
                             _NAD83_, _NTF_, _SPACE_, _WGS84_, _0_0, _0_5, _1_0, _90_0
from pygeodesy.interns import _C_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _lazyNamedEnumItem as _lazy, _NamedBase, \
                            _NamedEnum, _NamedEnumItem, nameof, _xnamed
from pygeodesy.namedTuples import EasNor3Tuple, LatLonDatum3Tuple, \
                                  LatLon2Tuple, _LL4Tuple, PhiLam2Tuple
from pygeodesy.props import deprecated_method, Property_RO, _update_all
from pygeodesy.streprs import Fmt, _fstrENH2, _xzipairs
from pygeodesy.units import Easting, Height, Lam_, Northing, Phi_, Scalar_
from pygeodesy.utily import degrees90, degrees180, sincos2, tanPI_2_2

from math import atan, log, radians, sin, sqrt

__all__ = _ALL_LAZY.lcc
__version__ = '21.12.30'

_E0_   = 'E0'
_N0_   = 'N0'
_par1_ = 'par1'
_par2_ = 'par2'
_SP_   = 'SP'


class Conic(_NamedEnumItem):
    '''Lambert conformal conic projection (1- or 2-SP).
    '''
    _auth  = NN    # authorization (C{str})
    _datum = None  # datum (L{Datum})
    _name  = NN    # conic (L{Conic})

    _e  = _0_0  # ellipsoid excentricity (C{float})
    _E0 = _0_0  # false easting (C{float})
    _k0 = _1_0  # scale factor (C{float})
    _N0 = _0_0  # false northing (C{float})
    _SP =  0    # 1- or 2-SP (C{int})

    _opt3 = _0_0  # optional, longitude (C{radians})
    _par1 = _0_0  # 1st std parallel (C{radians})
    _par2 = _0_0  # 2nd std parallel (C{radians})
    _phi0 = _0_0  # origin lat (C{radians})
    _lam0 = _0_0  # origin lon (C{radians})

    _aF = _0_0  # precomputed F (C{float})
    _n  = _0_0  # precomputed n (C{float})
    _n_ = _0_0  # precomputed 1 / n (C{float})
    _r0 = _0_0  # precomputed rho0 (C{float})

    def __init__(self, latlon0, par1, par2=None, E0=0, N0=0,
                       k0=1, opt3=0, name=NN, auth=NN):
        '''New Lambert conformal conic projection.

           @arg latlon0: Origin with (ellipsoidal) datum (C{LatLon}).
           @arg par1: First standard parallel (C{degrees90}).
           @kwarg par2: Optional, second standard parallel (C{degrees90}).
           @kwarg E0: Optional, false easting (C{meter}).
           @kwarg N0: Optional, false northing (C{meter}).
           @kwarg k0: Optional scale factor (C{scalar}).
           @kwarg opt3: Optional meridian (C{degrees180}).
           @kwarg name: Optional name of the conic (C{str}).
           @kwarg auth: Optional authentication authority (C{str}).

           @return: A Lambert projection (L{Conic}).

           @raise TypeError: Non-ellipsoidal B{C{latlon0}}.

           @raise ValueError: Invalid B{C{par1}}, B{C{par2}},
                              B{C{E0}}, B{C{N0}}, B{C{k0}}
                              or B{C{opt3}}.

           @example:

            >>> from pygeodesy import Conic, Datums, ellipsoidalNvector
            >>> ll0 = ellipsoidalNvector.LatLon(23, -96, datum=Datums.NAD27)
            >>> Snyder = Conic(ll0, 33, 45, E0=0, N0=0, name='Snyder')
        '''
        if latlon0 is not None:
            _xinstanceof(_LLEB, latlon0=latlon0)
            self._phi0, self._lam0 = latlon0.philam

            self._par1 = Phi_(par1=par1)
            self._par2 = self._par1 if par2 is None else Phi_(par2=par2)

            if k0 != 1:
                self._k0 = Scalar_(k0=k0)
            if E0:
                self._E0 = Northing(E0=E0, falsed=True)
            if N0:
                self._N0 = Easting(N0=N0, falsed=True)
            if opt3:
                self._opt3 = Lam_(opt3=opt3)

            self.toDatum(latlon0.datum)._dup2(self)
            self._register(Conics, name)
        elif name:
            self.name = name
        if auth:
            self._auth = str(auth)

    @Property_RO
    def auth(self):
        '''Get the authentication authority (C{str}).
        '''
        return self._auth

    @deprecated_method
    def convertDatum(self, datum):
        '''DEPRECATED, use method L{Conic.toDatum}.'''
        return self.toDatum(datum)

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def E0(self):
        '''Get the false easting (C{meter}).
        '''
        return self._E0

    @Property_RO
    def k0(self):
        '''Get scale factor (C{float}).
        '''
        return self._k0

    @Property_RO
    def lat0(self):
        '''Get the origin latitude (C{degrees90}).
        '''
        return degrees90(self._phi0)

    @Property_RO
    def latlon0(self):
        '''Get the central origin (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat0, self.lon0, name=self.name)

    @Property_RO
    def lam0(self):
        '''Get the central meridian (C{radians}).
        '''
        return self._lam0

    @Property_RO
    def lon0(self):
        '''Get the central meridian (C{degrees180}).
        '''
        return degrees180(self._lam0)

    @Property_RO
    def N0(self):
        '''Get the false northing (C{meter}).
        '''
        return self._N0

    @Property_RO
    def name2(self):
        '''Get the conic and datum names as "conic.datum" (C{str}).
        '''
        return self._DOT_(self.datum.name)

    @Property_RO
    def opt3(self):
        '''Get the optional meridian (C{degrees180}).
        '''
        return degrees180(self._opt3)

    @Property_RO
    def par1(self):
        '''Get the 1st standard parallel (C{degrees90}).
        '''
        return degrees90(self._par1)

    @Property_RO
    def par2(self):
        '''Get the 2nd standard parallel (C{degrees90}).
        '''
        return degrees90(self._par2)

    @Property_RO
    def phi0(self):
        '''Get the origin latitude (C{radians}).
        '''
        return self._phi0

    @Property_RO
    def philam0(self):
        '''Get the central origin (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(self.phi0, self.lam0, name=self.name)

    @Property_RO
    def SP(self):
        '''Get the number of standard parallels (C{int}).
        '''
        return self._SP

    def toDatum(self, datum):
        '''Convert this conic to the given datum.

           @arg datum: Ellipsoidal datum to use (L{Datum}, L{Ellipsoid},
                       L{Ellipsoid2} or L{a_f2Tuple}).

           @return: Converted conic, unregistered (L{Conic}).

           @raise TypeError: Non-ellipsoidal B{C{datum}}.
        '''
        d = _ellipsoidal_datum(datum, name=self.name)
        E =  d.ellipsoid
        if not E.isEllipsoidal:
            raise _IsnotError(_ellipsoidal_, datum=datum)

        c = self
        if c._e != E.e or c._datum != d:

            c = Conic(None, 0, name=self._name)
            self._dup2(c)
            c._datum = d
            c._e = E.e

            if abs(c._par1 - c._par2) < EPS:
                m1 = c._mdef(c._phi0)
                t1 = c._tdef(c._phi0)
                t0 = t1
                k  = 1  # _1_0
                n  = sin(c._phi0)
                sp = 1
            else:
                m1 = c._mdef(c._par1)
                m2 = c._mdef(c._par2)
                t1 = c._tdef(c._par1)
                t2 = c._tdef(c._par2)
                t0 = c._tdef(c._phi0)
                k  = c._k0
                n  = (log(m1) - log(m2)) \
                   / (log(t1) - log(t2))
                sp = 2

            F = m1 / (n * pow(t1, n))

            c._aF = k * E.a * F
            c._n  = n
            c._n_ = _1_0 / n
            c._r0 = c._rdef(t0)
            c._SP = sp

        return c

    def toStr(self, prec=8):  # PYCHOK expected
        '''Return this conic as a string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Conic attributes (C{str}).
        '''
        if self._SP == 1:
            return self._instr(prec, _lat0_, _lon0_, _par1_,
                                     _E0_, _N0_, _k0_, _SP_,
                                      datum=self.datum)
        else:
            return self._instr(prec, _lat0_, _lon0_, _par1_, _par2_,
                                     _E0_, _N0_, _k0_, _SP_,
                                      datum=self.datum)

    def _dup2(self, c):
        '''(INTERNAL) Copy this conic to C{c}.

           @arg c: Duplicate (L{Conic}).
        '''
        _update_all(c)

        c._auth  = self._auth
        c._datum = self._datum

        c._e  = self._e
        c._E0 = self._E0
        c._k0 = self._k0
        c._N0 = self._N0
        c._SP = self._SP

        c._par1 = self._par1
        c._par2 = self._par2
        c._phi0 = self._phi0
        c._lam0 = self._lam0
        c._opt3 = self._opt3

        c._aF = self._aF
        c._n  = self._n
        c._n_ = self._n_
        c._r0 = self._r0

    def _mdef(self, a):
        '''(INTERNAL) Compute m(a).
        '''
        s, c = sincos2(a)
        s = _1_0 - (s * self._e)**2
        return (c / sqrt(s)) if s > EPS02 else _0_0

    def _pdef(self, a):
        '''(INTERNAL) Compute p(a).
        '''
        s = self._e * sin(a)
        return pow((_1_0 - s) / (_1_0 + s), self._e * _0_5)

    def _rdef(self, t):
        '''(INTERNAL) Compute r(t).
        '''
        return self._aF * pow(t, self._n)

    def _tdef(self, a):
        '''(INTERNAL) Compute t(lat).
        '''
        return max(_0_0, tanPI_2_2(-a) / self._pdef(a))

    def _xdef(self, t_x):
        '''(INTERNAL) Compute x(t_x).
        '''
        return PI_2 - 2 * atan(t_x)  # XXX + self._phi0


Conic._name = Conic.__name__


class Conics(_NamedEnum):
    '''(INTERNAL) L{Conic} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, lat, lon, datum_name, *args, **kwds):
        '''(INTERNAL) Instantiate the L{Conic}.
        '''
        return Conic(_LLEB(lat, lon, datum=Datums.get(datum_name)), *args, **kwds)

Conics = Conics(Conic)  # PYCHOK singleton
'''Some pre-defined L{Conic}s, all I{lazily} instantiated.'''
Conics._assert(  # <https://SpatialReference.org/ref/sr-org/...>
#   AsLb   = _lazy('AsLb',   _F(-14.2666667), _F(170), _NAD27_, _0_0, _0_0,
#                             E0=_F(500000), N0=_0_0, auth='EPSG:2155'),  # American Samoa ... SP=1 !
    Be08Lb = _lazy('Be08Lb', _F(50.7978150), _F(4.359215833), _GRS80_, _F(49.8333339), _F(51.1666672),
                              E0=_F(649328.0), N0=_F(665262.0), auth='EPSG:9802'),  # Belgium
    Be72Lb = _lazy('Be72Lb', _90_0, _F(4.3674867), _NAD83_, _F(49.8333339), _F(51.1666672),
                              E0=_F(150000.013), N0=_F(5400088.438), auth='EPSG:31370'),  # Belgium
    Fr93Lb = _lazy('Fr93Lb', _F(46.5), _F(3), _WGS84_, _F(49), _F(44),
                              E0=_F(700000), N0=_F(6600000), auth='EPSG:2154'),  # RFG93, France
    MaNLb  = _lazy('MaNLb',  _F(33.3), _F(-5.4), _NTF_, _F(31.73), _F(34.87),
                              E0=_F(500000), N0=_F(300000)),  # Marocco
    MxLb   = _lazy('MxLb',   _F(12), _F(-102), _WGS84_, _F(17.5), _F(29.5),
                              E0=_F(2500000), N0=_0_0, auth='EPSG:2155'),  # Mexico
    PyT_Lb = _lazy('PyT_Lb', _F(46.8), _F(2.33722917), _NTF_, _F(45.89893890000052), _F(47.69601440000037),
                              E0=_F(600000), N0=_F(200000), auth='Test'),  # France?
    USA_Lb = _lazy('USA_Lb', _F(23), _F(-96), _WGS84_, _F(33), _F(45),
                              E0=_0_0, N0=_0_0),  # Conterminous, contiguous USA?
    WRF_Lb = _lazy('WRF_Lb', _F(40), _F(-97), _WGS84_, _F(33), _F(45),
                              E0=_0_0, N0=_0_0, auth='EPSG:4326')  # World
)


class LCCError(_ValueError):
    '''Lambert Conformal Conic C{LCC} or other L{Lcc} issue.
    '''
    pass


class Lcc(_NamedBase):
    '''Lambert conformal conic East-/Northing location.
    '''
    _conic    =  None  # Lambert projection (L{Conic})
    _easting  = _0_0   # Easting (C{float})
    _height   =  0     # height (C{meter})
    _northing = _0_0   # Northing (C{float})

    def __init__(self, e, n, h=0, conic=Conics.WRF_Lb, name=NN):
        '''New L{Lcc} Lamber conformal conic position.

           @arg e: Easting (C{meter}).
           @arg n: Northing (C{meter}).
           @kwarg h: Optional height (C{meter}).
           @kwarg conic: Optional, the conic projection (L{Conic}).
           @kwarg name: Optional name (C{str}).

           @return: The Lambert location (L{Lcc}).

           @raise LCCError: Invalid B{C{h}} or invalid or
                            negative B{C{e}} or B{C{n}}.

           @raise TypeError: If B{C{conic}} is not L{Conic}.

           @example:

            >>> lb = Lcc(448251, 5411932.0001)
        '''
        _xinstanceof(Conic, conic=conic)
        self._conic = conic
        self._easting  = Easting(e,  falsed=conic.E0 > 0, Error=LCCError)
        self._northing = Northing(n, falsed=conic.N0 > 0, Error=LCCError)
        if h:
            self._height = Height(h=h, Error=LCCError)
        if name:
            self.name = name

    @Property_RO
    def conic(self):
        '''Get the conic projection (L{Conic}).
        '''
        return self._conic

    @Property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @Property_RO
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude in C{degrees} (L{LatLon2Tuple}).
        '''
        ll = self.toLatLon(LatLon=None, datum=None)
        return LatLon2Tuple(ll.lat, ll.lon, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get the lat-, longitude and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

    @Property_RO
    def latlonheightdatum(self):
        '''Get the lat-, longitude in C{degrees} with height and datum (L{LatLon4Tuple}C{(lat, lon, height, datum)}).
        '''
        return self.latlonheight.to4Tuple(self.conic.datum)

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude in C{radians} (L{PhiLam2Tuple}).
        '''
        return PhiLam2Tuple(radians(self.latlon.lat),
                            radians(self.latlon.lon), name=self.name)

    @Property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    @Property_RO
    def philamheightdatum(self):
        '''Get the lat-, longitude in C{radians} with height and datum (L{PhiLam4Tuple}C{(phi, lam, height, datum)}).
        '''
        return self.philamheight.to4Tuple(self.datum)

    @deprecated_method
    def to3lld(self, datum=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{toLatLon}.

           @kwarg datum: Optional datum to use, otherwise use this
                         B{C{Lcc}}'s conic.datum (C{Datum}).

           @return: A L{LatLonDatum3Tuple}C{(lat, lon, datum)}.

           @raise TypeError: If B{C{datum}} is not ellipsoidal.
        '''
        if datum in (None, self.conic.datum):
            r = LatLonDatum3Tuple(self.latlon.lat,
                                  self.latlon.lon,
                                  self.conic.datum, name=self.name)
        else:
            r = self.toLatLon(LatLon=None, datum=datum)
            r = LatLonDatum3Tuple(r.lat, r.lon, r.datum, name=r.name)
        return r

    def toLatLon(self, LatLon=None, datum=None, height=None, **LatLon_kwds):
        '''Convert this L{Lcc} to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg datum: Optional datum to use, otherwise use this
                         B{C{Lcc}}'s conic.datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg height: Optional height for the point, overriding
                          the default height (C{meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The point (B{C{LatLon}}) or a
                    L{LatLon4Tuple}C{(lat, lon, height, datum)}
                    if B{C{LatLon}} is C{None}.

           @raise TypeError: If B{C{LatLon}} or B{C{datum}} is
                             not ellipsoidal or not valid.
        '''
        if LatLon:
            _xsubclassof(_LLEB, LatLon=LatLon)

        c = self.conic
        if datum not in (None, self.conic.datum):
            c = c.toDatum(datum)

        e =         self.easting  - c._E0
        n = c._r0 - self.northing + c._N0

        r_ = copysign0(hypot(e, n), c._n)
        t_ = pow(r_ / c._aF, c._n_)

        x = c._xdef(t_)  # XXX c._lam0
        while True:
            p, x = x, c._xdef(t_ * c._pdef(x))
            if abs(x - p) < 1e-9:  # XXX EPS too small?
                break
        lat = degrees90(x)
        lon = degrees180((atan(e / n) + c._opt3) * c._n_ + c._lam0)

        h = self.height if height is None else Height(height)
        return _LL4Tuple(lat, lon, h, c.datum, LatLon, LatLon_kwds,
                                               inst=self, name=self.name)

    def toRepr(self, prec=0, fmt=Fmt.SQUARE, sep=_COMMASPACE_, m=_m_, C=False, **unused):  # PYCHOK expected
        '''Return a string representation of this L{Lcc} position.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional, enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:values (C{str}).
           @kwarg m: Optional unit of the height, default meter (C{str}).
           @kwarg C: Optionally, include name of conic and datum (C{bool}).

           @return: This Lcc as "[E:meter, N:meter, H:m, C:Conic.Datum]"
                   (C{str}).
        '''
        t, T = _fstrENH2(self, prec, m)
        if C:
            t +=  self.conic.name2,
            T += _C_,
        return _xzipairs(T, t, sep=sep, fmt=fmt)

    def toStr(self, prec=0, sep=_SPACE_, m=_m_):  # PYCHOK expected
        '''Return a string representation of this L{Lcc} position.

           @kwarg prec: Optional number of decimal, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg m: Optional height units, default C{meter} (C{str}).

           @return: This Lcc as I{"easting nothing"} in C{meter} plus
                    I{" height"} suffixed with B{C{m}} if height is
                    non-zero (C{str}).
        '''
        t, _ = _fstrENH2(self, prec, m)
        return t if sep is None else sep.join(t)


def toLcc(latlon, conic=Conics.WRF_Lb, height=None, Lcc=Lcc, name=NN,
                                                  **Lcc_kwds):
    '''Convert an (ellipsoidal) geodetic point to a I{Lambert} location.

       @arg latlon: Ellipsoidal point (C{LatLon}).
       @kwarg conic: Optional Lambert projection to use (L{Conic}).
       @kwarg height: Optional height for the point, overriding the
                      default height (C{meter}).
       @kwarg Lcc: Optional class to return the I{Lambert} location
                   (L{Lcc}).
       @kwarg name: Optional B{C{Lcc}} name (C{str}).
       @kwarg Lcc_kwds: Optional, additional B{C{Lcc}} keyword
                        arguments, ignored if B{C{Lcc}} is C{None}.

       @return: The I{Lambert} location (L{Lcc}) or an
                L{EasNor3Tuple}C{(easting, northing, height)}
                if C{B{Lcc} is None}.

       @raise TypeError: If B{C{latlon}} is not ellipsoidal.
    '''
    _xinstanceof(_LLEB, latlon=latlon)

    a, b = latlon.philam
    c = conic.toDatum(latlon.datum)

    t = c._n * (b - c._lam0) - c._opt3
    st, ct = sincos2(t)

    r = c._rdef(c._tdef(a))
    e = c._E0         + r * st
    n = c._N0 + c._r0 - r * ct

    h = latlon.height if height is None else Height(height)
    r = EasNor3Tuple(e, n, h) if Lcc is None else \
                 Lcc(e, n, h=h, conic=c, **Lcc_kwds)
    return _xnamed(r, name or nameof(latlon))


if __name__ == '__main__':

    from pygeodesy.interns import _NL_, _NL_var_

    # __doc__ of this file, force all into registery
    t = _NL_ + Conics.toRepr(all=True)
    print(_NL_var_.join(t.split(_NL_)))

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
