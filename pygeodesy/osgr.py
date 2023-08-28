
# -*- coding: utf-8 -*-

u'''Ordnance Survey Grid References (OSGR) references on the UK U{National Grid
<https://www.OrdnanceSurvey.co.UK/documents/resources/guide-to-nationalgrid.pdf>}.

Classes L{Osgr} and L{OSGRError} and functions L{parseOSGR} and L{toOsgr}.

A pure Python implementation, transcoded from I{Chris Veness}' JavaScript originals U{OS National Grid
<https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>} and U{Module osgridref
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-osgridref.html>} and I{Charles Karney}'s
C++ class U{OSGB<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1OSGB.html>}.

OSGR provides geocoordinate references for UK mapping purposes, converted in 2015 to work with the C{WGS84}
or the original C{OSGB36} datum.  In addition, this implementation includes both the OS recommended and the
Krüger-based method to convert between OSGR and geodetic coordinates (with keyword argument C{kTM} of
function L{toOsgr}, method L{Osgr.toLatLon} and method C{toOsgr} of any ellipsoidal C{LatLon} class).

See U{Transverse Mercator: Redfearn series<https://WikiPedia.org/wiki/Transverse_Mercator:_Redfearn_series>},
Karney's U{"Transverse Mercator with an accuracy of a few nanometers", 2011<https://ArXiv.org/pdf/1002.1417v3.pdf>}
(building on U{"Konforme Abbildung des Erdellipsoids in der Ebene", 1912<https://bib.GFZ-Potsdam.DE/pub/digi/krueger2.pdf>},
U{"Die Mathematik der Gauß-Krueger-Abbildung", 2006<https://DE.WikiPedia.org/wiki/Gauß-Krüger-Koordinatensystem>},
U{A Guide<https://www.OrdnanceSurvey.co.UK/documents/resources/guide-coordinate-systems-great-britain.pdf>}
and U{Ordnance Survey National Grid<https://WikiPedia.org/wiki/Ordnance_Survey_National_Grid>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import halfs2, isbool, isfloat, map1, \
                            _splituple, _xsubclassof
from pygeodesy.constants import _1_0, _10_0,  _N_2_0  # PYCHOK used!
from pygeodesy.datums import Datums, _ellipsoidal_datum, _WGS84
# from pygeodesy.dms import parseDMS2   # _MODS
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _parseX, _TypeError, _ValueError, \
                             _xkwds, _xkwds_get
from pygeodesy.fmath import Fdot, fpowers, Fsum
# from pygeodesy.fsums import Fsum  # from .fmath
from pygeodesy.interns import MISSING, NN, _A_, _COLON_, _COMMA_, \
                             _COMMASPACE_, _DOT_, _ellipsoidal_, \
                             _latlon_, _not_, _SPACE_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, nameof, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, LatLon2Tuple, \
                                  LatLonDatum3Tuple
from pygeodesy.props import Property_RO, property_RO
from pygeodesy.streprs import _EN_WIDE, enstr2, _enstr2m3, Fmt, \
                              _resolution10, unstr, _xzipairs
from pygeodesy.units import Easting, Lam_, Lat, Lon, Northing, \
                            Phi_, Scalar, _10um, _100km
from pygeodesy.utily import degrees90, degrees180, sincostan3, truncate

from math import cos, fabs, radians, sin, sqrt

__all__ = _ALL_LAZY.osgr
__version__ = '23.08.24'

_equivalent_ = 'equivalent'
_OSGR_       = 'OSGR'
_ord_A       =  ord(_A_)
_TRIPS       =  33  # .toLatLon convergence


class _NG(object):
    '''Ordnance Survey National Grid parameters.
    '''
    @Property_RO
    def a0(self):  # equatoradius, scaled
        return self.ellipsoid.a * self.k0

    @Property_RO
    def b0(self):  # polaradius, scaled
        return self.ellipsoid.b * self.k0

    @Property_RO
    def datum(self):  # datum, Airy130 ellipsoid
        return Datums.OSGB36

    @Property_RO
    def eas0(self):  # False origin easting (C{meter})
        return Easting(4 * _100km)

    @Property_RO
    def easX(self):  # easting [0..extent] (C{meter})
        return Easting(7 * _100km)

    @Property_RO
    def ellipsoid(self):  # ellipsoid, Airy130
        return self.datum.ellipsoid

    def forward2(self, latlon):  # convert C{latlon} to (easting, norting), as I{Karney}'s
        # U{Forward<https://GeographicLib.SourceForge.io/C++/doc/OSGB_8hpp_source.html>}
        t = self.kTM.forward(latlon.lat, latlon.lon, lon0=self.lon0)
        e = t.easting  + self.eas0
        n = t.northing + self.nor0ffset
        return e, n

    @Property_RO
    def k0(self):  # central scale (C{float}), like I{Karney}'s CentralScale
        # <https://GeographicLib.SourceForge.io/C++/doc/OSGB_8hpp_source.html>
        _0_9998268 = (9998268 - 10000000) / 10000000
        return Scalar(_10_0**_0_9998268)  # 0.9996012717...

    @Property_RO
    def kTM(self):  # the L{KTransverseMercator} instance, like I{Karney}'s OSGBTM
        # <https://GeographicLib.SourceForge.io/C++/doc/OSGB_8cpp_source.html>
        return _MODS.ktm.KTransverseMercator(self.datum, lon0=0, k0=self.k0)

    @Property_RO
    def lam0(self):  # True origin longitude C{radians}
        return Lam_(self.lon0)

    @Property_RO
    def lat0(self):  # True origin latitude, 49°N
        return Lat(49.0)

    @Property_RO
    def lon0(self):  # True origin longitude, 2°W
        return Lon(_N_2_0)

    @Property_RO
    def Mabcd(self):  # meridional coefficients (a, b, c, d)
        n, n2, n3 = fpowers(self.ellipsoid.n, 3)
        M = (Fsum(4,  4 * n,  5 * n2,  5 * n3) / 4,
             Fsum(   24 * n, 24 * n2, 21 * n3) / 8,
             Fsum(           15 * n2, 15 * n3) / 8,
                                     (35 * n3 / 24))
        return M

    def Mabcd0(self, a):  # meridional arc, scaled
        c = a + self.phi0
        s = a - self.phi0
        R = Fdot(self.Mabcd, s, -sin(s)     * cos(c),
                                 sin(s * 2) * cos(c * 2),
                                -sin(s * 3) * cos(c * 3))
        return float(R * self.b0)

    @Property_RO
    def nor0(self):  # False origin northing (C{meter})
        return Northing(-_100km)

    @Property_RO
    def nor0ffset(self):  # like I{Karney}'s computenorthoffset
        # <https://GeographicLib.SourceForge.io/C++/doc/OSGB_8cpp_source.html>
        return self.nor0 - self.kTM.forward(self.lat0, 0).northing

    @Property_RO
    def norX(self):  # northing [0..extent] (C{meter})
        return Northing(13 * _100km)

    def nu_rho_eta3(self, sa):  # 3-tuple (nu, nu / rho, eta2)
        E = self.ellipsoid  # rho, nu = E.roc2_(sa)  # .k0?
        s = E.e2s2(sa)  # == 1 - E.e2 * sa**2
        v = self.a0 / sqrt(s)  # == nu, transverse roc
        # rho = .a0 * E.e21 / s**1.5 == v * E.e21 / s
        # r = v * E.e21 / s  # == rho, meridional roc
        # nu / rho == v / (v * E.e21 / s) == s / E.e21 == ...
        s *= E._1_e21  # ... s * E._1_e21 == s * E.a2_b2
        return v, s, (s - _1_0)  # η2 = nu / rho - 1

    @Property_RO
    def phi0(self):  # True origin latitude C{radians}
        return Phi_(self.lat0)

    def reverse(self, osgr):  # convert C{osgr} to (ellipsoidal} LatLon, as I{Karney}'s
        # U{Reverse<https://GeographicLib.SourceForge.io/C++/doc/OSGB_8hpp_source.html>}
        r = osgr._latlonTM
        if r is None:
            x =  osgr.easting  - self.eas0
            y =  osgr.northing - self.nor0ffset
            t =  self.kTM.reverse(x, y, lon0=self.lon0)
            r = _LLEB(t.lat, t.lon, datum=self.datum, name=osgr.name)
            osgr._latlonTM = r
        return r

_NG = _NG()  # PYCHOK singleton


class OSGRError(_ValueError):
    '''Error raised for a L{parseOSGR}, L{Osgr} or other OSGR issue.
    '''
    pass


class Osgr(_NamedBase):
    '''Ordnance Survey Grid References (OSGR) coordinates on
       the U{National Grid<https://www.OrdnanceSurvey.co.UK/
       documents/resources/guide-to-nationalgrid.pdf>}.
    '''
    _datum      = _NG.datum  # default datum (L{Datums.OSGB36})
    _easting    =  0         # Easting (C{meter})
    _latlon     =  None      # cached B{C{_toLatlon}}
    _latlonTM   =  None      # cached B{C{_toLatlon kTM}}
    _northing   =  0         # Nothing (C{meter})
    _resolution =  0         # from L{parseOSGR} (C{meter})

    def __init__(self, easting, northing, datum=None, name=NN,
                                          resolution=0):
        '''New L{Osgr} coordinate.

           @arg easting: Easting from the OS C{National Grid}
                         origin (C{meter}).
           @arg northing: Northing from the OS C{National Grid}
                          origin (C{meter}).
           @kwarg datum: Override default datum (C{Datums.OSGB36}).
           @kwarg name: Optional name (C{str}).
           @kwarg resolution: Optional resolution (C{meter}),
                              C{0} for default.

           @raise OSGRError: Invalid or negative B{C{easting}} or
                             B{C{northing}} or B{C{datum}} not an
                             C{Datums.OSGB36} equivalent.
        '''
        if datum:  # PYCHOK no cover
            try:
                self._datum = _ellipsoidal_datum(datum)
                if self.datum != _NG.datum:
                    raise ValueError(_not_(_NG.datum.name, _equivalent_))
            except (TypeError, ValueError) as x:
                raise OSGRError(datum=datum, cause=x)

        self._easting  = Easting( easting,  Error=OSGRError, high=_NG.easX)
        self._northing = Northing(northing, Error=OSGRError, high=_NG.norX)

        if name:
            self.name = name
        if resolution:
            self._resolution = _resolution10(resolution, Error=OSGRError)

    def __str__(self):
        return self.toStr(GD=True, sep=_SPACE_)

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @Property_RO
    def falsing0(self):
        '''Get the C{OS National Grid} falsing (L{EasNor2Tuple}).
        '''
        return EasNor2Tuple(_NG.eas0, _NG.nor0, name=_OSGR_)

    @property_RO
    def iteration(self):
        '''Get the most recent C{Osgr.toLatLon} iteration number
           (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @Property_RO
    def latlon0(self):
        '''Get the C{OS National Grid} origin (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(_NG.lat, _NG.lon0, name=_OSGR_)

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    def parse(self, strOSGR, name=NN):
        '''Parse an OSGR reference to a similar L{Osgr} instance.

           @arg strOSGR: The OSGR reference (C{str}), see function L{parseOSGR}.
           @kwarg name: Optional instance name (C{str}), overriding this name.

           @return: The similar instance (L{Osgr})

           @raise OSGRError: Invalid B{C{strOSGR}}.
        '''
        return parseOSGR(strOSGR, Osgr=self.classof, name=name or self.name)

    @property_RO
    def resolution(self):
        '''Get the OSGR resolution (C{meter}, power of 10) or C{0} if undefined.
        '''
        return self._resolution

    def toLatLon(self, LatLon=None, datum=_WGS84, kTM=False, eps=_10um, **LatLon_kwds):
        '''Convert this L{Osgr} coordinate to an (ellipsoidal) geodetic
           point.

           @kwarg LatLon: Optional ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg datum: Optional datum to convert to (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{Ellipsoid2}
                         or L{a_f2Tuple}).
           @kwarg kTM: If C{True} use I{Karney}'s Krüger method from
                       module L{ktm}, otherwise use the Ordnance
                       Survey formulation (C{bool}).
           @kwarg eps: Tolerance for OS convergence (C{meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: A B{C{LatLon}} instance or if B{C{LatLon}} is C{None}
                    a L{LatLonDatum3Tuple}C{(lat, lon, datum)}.

           @note: While OS grid references are based on the OSGB36 datum,
                  the Ordnance Survey have deprecated the use of OSGB36 for
                  lat-/longitude coordinates (in favour of WGS84).  Hence,
                  this method returns WGS84 by default with OSGB36 as an
                  option, U{see<https://www.OrdnanceSurvey.co.UK/blog/2014/12/2>}.

           @note: The formulation implemented here due to Thomas, Redfearn,
                  etc. is as published by the Ordnance Survey, but is
                  inferior to Krüger as used by e.g. Karney 2011.

           @raise OSGRError: No convergence.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal or B{C{datum}}
                             is invalid or conversion to B{C{datum}} failed.

           @example:

            >>> from pygeodesy import Datums, ellipsoidalVincenty as eV, Osgr
            >>> g = Osgr(651409.903, 313177.270)
            >>> p = g.toLatLon(eV.LatLon)  # 52°39′28.723″N, 001°42′57.787″E
            >>> p = g.toLatLon(eV.LatLon, kTM=True)  # 52°39′28.723″N, 001°42′57.787″E
            >>> # to obtain (historical) OSGB36 lat-/longitude point
            >>> p = g.toLatLon(eV.LatLon, datum=Datums.OSGB36)  # 52°39′27.253″N, 001°43′04.518″E
        '''
        NG = _NG
        if kTM:
            r = NG.reverse(self)

        elif self._latlon is None:
            e0 =     self.easting  - NG.eas0
            n0 = m = self.northing - NG.nor0

            _M = NG.Mabcd0
            a0 = NG.a0
            a  = NG.phi0
            _A = Fsum(a).fsum_
            for self._iteration in range(1, _TRIPS):
                a = _A(m / a0)
                m = n0 - _M(a)  # meridional arc
                if fabs(m) < eps:
                    break
            else:  # PYCHOK no cover
                t =  str(self)
                t =  Fmt.PAREN(self.classname, repr(t))
                t = _DOT_(t, self.toLatLon.__name__)
                t =  unstr(t, eps=eps, kTM=kTM)
                raise OSGRError(Fmt.no_convergence(m), txt=t)

            sa, ca, ta = sincostan3(a)
            v, v_r, n2 = NG.nu_rho_eta3(sa)

            ta2 = ta**2
            ta4 = ta2**2

            ta *= v_r / 2
            d   = e0 / v
            d2  = d**2

            a = (d2 * ta * (-1 +  # Horner-like
                 d2 / 12 * (Fsum( 5,  3 * ta2, -9 * ta2 * n2, n2) -
                 d2 / 30 *  Fsum(61, 90 * ta2, 45 * ta4)))).fsum_(a)

            b = (d  / ca * (1 -  # Horner-like
                 d2 /  6 * (Fsum(v_r,  2 * ta2) -
                 d2 / 20 * (Fsum( 5,  28 * ta2,   24 * ta4) +
                 d2 / 42 *  Fsum(61, 662 * ta2, 1320 * ta4,
                                     720 * ta2 * ta4))))).fsum_(NG.lam0)

            r = _LLEB(degrees90(a), degrees180(b), datum=self.datum, name=self.name)
            r._iteration = self._iteration  # only ellipsoidal LatLon
            self._latlon = r
        else:
            r = self._latlon

        return _ll2LatLon3(r, LatLon, datum, LatLon_kwds)

    @Property_RO
    def scale0(self):
        '''Get the C{OS National Grid} central scale (C{scalar}).
        '''
        return _NG.k0

    def toRepr(self, GD=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **prec):  # PYCHOK expected
        '''Return a string representation of this L{Osgr} coordinate.

           @kwarg GD: If C{bool}, in- or exclude the 2-letter grid designation and get
                      the new B{C{prec}} behavior, otherwise if C{None}, default to the
                      DEPRECATED definition C{B{prec}=5} I{for backward compatibility}.
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Separator to join (C{str}) or C{None} to return an unjoined 2- or
                       3-C{tuple} of C{str}s.
           @kwarg prec: Precison C{B{prec}=0}, the number of I{decimal} digits (C{int}) or
                        if negative, the number of I{units to drop}, like MGRS U{PRECISION
                        <https://GeographicLib.SourceForge.io/C++/doc/GeoConvert.1.html#PRECISION>}.

           @return: This OSGR as (C{str}), C{"[G:GD, E:meter, N:meter]"} or if C{B{GD}=False}
                    C{"[OSGR:easting,northing]"} or C{B{GD}=False} and C{B{prec} > 0} if
                    C{"[OSGR:easting.d,northing.d]"}.

           @note: OSGR easting and northing values are truncated, not rounded.

           @raise OSGRError: If C{B{GD} not in (None, True, False)} or if C{B{prec} < -4}
                             and C{B{GD}=False}.

           @raise ValueError: Invalid B{C{prec}}.
        '''
        GD, prec = _GD_prec2(GD, fmt=fmt, sep=sep, **prec)

        if GD:
            t =  self.toStr(GD=True, prec=prec, sep=None)
            t = _xzipairs('GEN', t, sep=sep, fmt=fmt)
        else:
            t = _COLON_(_OSGR_, self.toStr(GD=False, prec=prec))
            if fmt:
                t = fmt % (t,)
        return t

    def toStr(self, GD=None, sep=NN, **prec):  # PYCHOK expected
        '''Return this L{Osgr} coordinate as a string.

           @kwarg GD: If C{bool}, in- or exclude the 2-letter grid designation and get
                      the new B{C{prec}} behavior, otherwise if C{None}, default to the
                      DEPRECATED definition C{B{prec}=5} I{for backward compatibility}.
           @kwarg sep: Separator to join (C{str}) or C{None} to return an unjoined 2- or
                       3-C{tuple} of C{str}s.
           @kwarg prec: Precison C{B{prec}=0}, the number of I{decimal} digits (C{int}) or
                        if negative, the number of I{units to drop}, like MGRS U{PRECISION
                        <https://GeographicLib.SourceForge.io/C++/doc/GeoConvert.1.html#PRECISION>}.

           @return: This OSGR as (C{str}), C{"GD meter meter"} or if C{B{GD}=False}
                    C{"easting,northing"} or if C{B{GD}=False} and C{B{prec} > 0}
                    C{"easting.d,northing.d"}

           @note: OSGR easting and northing values are truncated, not rounded.

           @raise OSGRError: If C{B{GD} not in (None, True, False)} or if C{B{prec}
                             < -4} and C{B{GD}=False}.

           @raise ValueError: Invalid B{C{prec}}.

           @example:

            >>> r = Osgr(651409, 313177)
            >>> str(r)  # 'TG 5140 1317'
            >>> r.toStr()  # 'TG5140913177'
            >>> r.toStr(GD=False)  # '651409,313177'
        '''
        def _i2c(i):
            if i > 7:
                i += 1
            return chr(_ord_A + i)

        GD, prec = _GD_prec2(GD, sep=sep, **prec)

        if GD:
            E, e = divmod(self.easting,  _100km)
            N, n = divmod(self.northing, _100km)
            E, N = int(E), int(N)
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                raise OSGRError(E=E, e=e, N=N, n=n, prec=prec, sep=sep)
            N  =  19 - N
            EN = _i2c( N - (N % 5) + (E + 10) // 5) + \
                 _i2c((N * 5) % 25 + (E % 5))
            t = enstr2(e, n, prec, EN)
            s = sep

        elif prec <= -_EN_WIDE:
            raise OSGRError(GD=GD, prec=prec, sep=sep)
        else:
            t = enstr2(self.easting, self.northing, prec, dot=True,
                                                    wide=_EN_WIDE + 1)
            s = sep if sep is None else (sep or _COMMA_)

        return t if s is None else s.join(t)


def _GD_prec2(GD, **prec_et_al):
    '''(INTERNAL) Handle C{prec} backward compatibility.
    '''
    if GD is None:  # old C{prec} 5+ or neg
        prec = _xkwds_get(prec_et_al, prec=_EN_WIDE)
        GD   =  prec > 0
        prec = (prec - _EN_WIDE) if GD else -prec
    elif isbool(GD):
        prec = _xkwds_get(prec_et_al, prec=0)
    else:
        raise OSGRError(GD=GD, **prec_et_al)
    return GD, prec


def _ll2datum(ll, datum, name):
    '''(INTERNAL) Convert datum if needed.
    '''
    if datum:
        try:
            if ll.datum != datum:
                ll = ll.toDatum(datum)
        except (AttributeError, TypeError, ValueError) as x:
            raise _TypeError(cause=x, datum=datum.name, **{name: ll})
    return ll


def _ll2LatLon3(ll, LatLon, datum, LatLon_kwds):
    '''(INTERNAL) Convert C{ll} to C{LatLon}
    '''
    if LatLon is None:
        r = _ll2datum(ll, datum, LatLonDatum3Tuple.__name__)
        r =  LatLonDatum3Tuple(r.lat, r.lon, r.datum)
    else:  # must be ellipsoidal
        _xsubclassof(_LLEB, LatLon=LatLon)
        r = _ll2datum(ll, datum, LatLon.__name__)
        r =  LatLon(r.lat, r.lon, datum=r.datum, **LatLon_kwds)
    if r._iteration != ll._iteration:
        r._iteration = ll._iteration
    return _xnamed(r, nameof(ll))


def parseOSGR(strOSGR, Osgr=Osgr, name=NN, **Osgr_kwds):
    '''Parse a string representing an OS Grid Reference, consisting
       of C{"[GD] easting northing"}.

       Accepts standard OS Grid References like "SU 387 148", with
       or without whitespace separators, from 2- up to 22-digit
       references, or all-numeric, comma-separated references in
       meters, for example "438700,114800".

       @arg strOSGR: An OSGR coordinate (C{str}).
       @kwarg Osgr: Optional class to return the OSGR coordinate
                    (L{Osgr}) or C{None}.
       @kwarg name: Optional B{C{Osgr}} name (C{str}).
       @kwarg Osgr_kwds: Optional, additional B{C{Osgr}} keyword
                         arguments, ignored if C{B{Osgr} is None}.

       @return: An (B{C{Osgr}}) instance or if B{C{Osgr}} is
                C{None} an L{EasNor2Tuple}C{(easting, northing)}.

       @raise OSGRError: Invalid B{C{strOSGR}}.

       @example:

        >>> r = parseOSGR('TG5140913177')
        >>> str(r)  # 'TG 51409 13177'
        >>> r = parseOSGR('TG 51409 13177')
        >>> r.toStr()  # 'TG5140913177'
        >>> r = parseOSGR('651409,313177')
        >>> r.toStr(sep=' ')  # 'TG 51409 13177'
        >>> r.toStr(GD=False)  # '651409,313177'
    '''
    def _c2i(G):
        g = ord(G.upper()) - _ord_A
        if g > 7:
            g -= 1
        if g < 0 or g > 25:
            raise ValueError
        return g

    def _OSGR(strOSGR, Osgr, name):
        s = _splituple(strOSGR.strip())
        p =  len(s)
        if not p:
            raise ValueError
        g = s[0]
        if p == 2 and isfloat(g):  # "easting,northing"
            e, n, m = _enstr2m3(*s, wide=_EN_WIDE + 1)

        else:
            if p == 1:  # "GReastingnorthing"
                s = halfs2(g[2:])
                g = g[:2]
            elif p == 2:  # "GReasting northing"
                s = g[2:], s[1]  # for backward ...
                g = g[:2]  # ... compatibility
            elif p != 3:
                raise ValueError
            else:  # "GR easting northing"
                s = s[1:]

            e, n = map(_c2i, g)
            n, m = divmod(n, 5)
            E = ((e - 2) % 5) * 5 + m
            N = 19 - (e // 5) * 5 - n
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                raise ValueError

            e, n, m = _enstr2m3(*s, wide=_EN_WIDE)
            e +=  E * _100km
            n +=  N * _100km

        if Osgr is None:
            _ = _MODS.osgr.Osgr(e, n, resolution=m)  # validate
            r =  EasNor2Tuple(e, n, name=name)
        else:
            r =  Osgr(e, n, name=name,
                         **_xkwds(Osgr_kwds, resolution=m))
        return r

    return _parseX(_OSGR, strOSGR, Osgr, name,
                          strOSGR=strOSGR, Error=OSGRError)


def toOsgr(latlon, lon=None, kTM=False, datum=_WGS84, Osgr=Osgr, name=NN,  # MCCABE 14
                                               **prec_Osgr_kwds):
    '''Convert a lat-/longitude point to an OSGR coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal) geodetic
                    C{LatLon} point.
       @kwarg lon: Optional longitude in degrees (scalar or C{None}).
       @kwarg kTM: If C{True} use I{Karney}'s Krüger method from
                   module L{ktm}, otherwise use the Ordnance
                   Survey formulation (C{bool}).
       @kwarg datum: Optional datum to convert B{C{lat, lon}} from
                     (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                     L{a_f2Tuple}).
       @kwarg Osgr: Optional class to return the OSGR coordinate
                    (L{Osgr}) or C{None}.
       @kwarg name: Optional B{C{Osgr}} name (C{str}).
       @kwarg prec_Osgr_kwds: Optional L{truncate} precision
                              C{B{prec}=ndigits} and/or additional
                              B{C{Osgr}} keyword arguments, ignored
                              if C{B{Osgr} is None}.

       @return: An (B{C{Osgr}}) instance or if B{C{Osgr}} is C{None}
                an L{EasNor2Tuple}C{(easting, northing)}.

       @note: If L{isint}C{(B{prec})} both easting and northing are
              L{truncate}d to the given number of digits.

       @raise OSGRError: Invalid B{C{latlon}} or B{C{lon}}.

       @raise TypeError: Non-ellipsoidal B{C{latlon}} or invalid
                         B{C{datum}}, B{C{Osgr}}, B{C{Osgr_kwds}}
                         or conversion to C{Datums.OSGB36} failed.

       @example:

        >>> p = LatLon(52.65798, 1.71605)
        >>> r = toOsgr(p)  # [G:TG, E:51409, N:13177]
        >>> # for conversion of (historical) OSGB36 lat-/longitude:
        >>> r = toOsgr(52.65757, 1.71791, datum=Datums.OSGB36)
        >>> # alternatively and using Krüger
        >>> r = p.toOsgr(kTM=True)  # [G:TG, E:51409, N:13177]
    '''
    def _prec_kwds2(prec=MISSING, **kwds):
        return prec, kwds

    if lon is not None:
        try:
            lat, lon = _MODS.dms.parseDMS2(latlon, lon)
            latlon   = _LLEB(lat, lon, datum=datum)
        except Exception as x:
            raise OSGRError(latlon=latlon, lon=lon, datum=datum, cause=x)
    elif not isinstance(latlon, _LLEB):
        raise _TypeError(latlon=latlon, txt=_not_(_ellipsoidal_))
    elif not name:  # use latlon.name
        name = nameof(latlon)

    NG = _NG
    # convert latlon to OSGB36 first
    ll = _ll2datum(latlon, NG.datum, _latlon_)

    if kTM:
        e, n = NG.forward2(ll)

    else:
        try:
            a, b = ll.philam
        except AttributeError:
            a, b = map1(radians, ll.lat, ll.lon)

        sa, ca, ta = sincostan3(a)
        v, v_r, n2 = NG.nu_rho_eta3(sa)

        m0  = NG.Mabcd0(a)
        b  -= NG.lam0
        t   = b * sa * v / 2
        d   = b * ca
        d2  = d**2

        ta2 = -(ta**2)
        ta4 =   ta2**2

        e = (d  *  v * (1 +  # Horner-like
             d2 /  6 * (Fsum(v_r, ta2) +
             d2 / 20 *  Fsum(5,  18 * ta2, ta4, 14 * n2,
                            58 * n2 * ta2)))).fsum_(NG.eas0)

        n = (d  *  t * (1 +  # Horner-like
             d2 / 12 * (Fsum( 5, ta2,  9 * n2) +
             d2 / 30 *  Fsum(61, ta4, 58 * ta2)))).fsum_(m0, NG.nor0)

    p, kwds = _prec_kwds2(**prec_Osgr_kwds)
    if p is not MISSING:
        e = truncate(e, p)
        n = truncate(n, p)

    if Osgr is None:
        _ = _MODS.osgr.Osgr(e, n)  # validate
        r =  EasNor2Tuple(e, n)
    else:
        r = Osgr(e, n, **kwds)  # datum=NG.datum
        if lon is None and isinstance(latlon, _LLEB):
            if kTM:
                r._latlonTM = latlon  # XXX weakref(latlon)?
            else:
                r._latlon = latlon  # XXX weakref(latlon)?
    return _xnamed(r, name or nameof(latlon))


if __name__ == '__main__':

    from pygeodesy.lazily import printf
    from random import random, seed
    from time import localtime

    seed(localtime().tm_yday)

    def _rnd(X, n):
        X -= 2
        d = set()
        while len(d) < n:
            r = 1 + int(random() * X)
            if r not in d:
                d.add(r)
                yield r

    D  = _NG.datum
    i  = t = 0
    t1 = t2 = 0, 0, 0, 0
    for e in _rnd(_NG.easX, 256):
        for n in _rnd(_NG.norX, 512):
            p  = False
            t += 1

            g = Osgr(e, n)
            v = g.toLatLon(kTM=False, datum=D)
            k = g.toLatLon(kTM=True,  datum=D)
            d = max(fabs(v.lat - k.lat), fabs(v.lon - k.lon))
            if d > t1[2]:
                t1 = e, n, d, t
                p  = True

            ll = _LLEB((v.lat + k.lat) / 2,
                       (v.lon + k.lon) / 2, datum=D)
            v  =  ll.toOsgr(kTM=False)
            k  =  ll.toOsgr(kTM=True)
            d  =  max(fabs(v.easting  - k.easting),
                      fabs(v.northing - k.northing))
            if d > t2[2]:
                t2 = ll.lat, ll.lon, d, t
                p  = True

            if p:
                i += 1
                printf('%5d: %s  %s', i,
                       'll(%.2f, %.2f) %.3e %d' % t2,
                       'en(%d, %d) %.3e %d' % t1)
    printf('%d total %s', t, D.name)

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

# % python3 -m pygeodesy.osgr
#     1: ll(53.42, -0.59) 4.672e-07 1  en(493496, 392519) 2.796e-11 1
#     2: ll(60.86, -0.28) 2.760e-05 2  en(493496, 1220986) 2.509e-10 2
#     3: ll(61.41, -0.25) 3.045e-05 13  en(493496, 1281644) 2.774e-10 13
#     4: ll(61.41, -0.25) 3.045e-05 13  en(493496, 1192797) 3.038e-10 20
#     5: ll(61.41, -0.25) 3.045e-05 13  en(493496, 1192249) 3.073e-10 120
#     6: ll(61.55, -0.24) 3.120e-05 160  en(493496, 1192249) 3.073e-10 120
#     7: ll(61.55, -0.24) 3.122e-05 435  en(493496, 1192249) 3.073e-10 120
#     8: ll(61.57, -0.24) 3.130e-05 473  en(493496, 1192249) 3.073e-10 120
#     9: ll(58.66, -8.56) 8.084e-04 513  en(19711, 993800) 3.020e-06 513
#    10: ll(52.83, -7.65) 8.156e-04 518  en(19711, 993800) 3.020e-06 513
#    11: ll(51.55, -7.49) 8.755e-04 519  en(19711, 993800) 3.020e-06 513
#    12: ll(60.20, -8.87) 9.439e-04 521  en(19711, 1165686) 4.318e-06 521
#    13: ll(60.45, -8.92) 9.668e-04 532  en(19711, 1194002) 4.588e-06 532
#    14: ll(61.17, -9.08) 1.371e-03 535  en(19711, 1274463) 5.465e-06 535
#    15: ll(61.31, -9.11) 1.463e-03 642  en(19711, 1290590) 5.663e-06 642
#    16: ll(61.35, -9.12) 1.488e-03 807  en(19711, 1294976) 5.718e-06 807
#    17: ll(61.38, -9.13) 1.510e-03 929  en(19711, 1298667) 5.765e-06 929
#    18: ll(61.11, -9.24) 1.584e-03 11270  en(10307, 1268759) 6.404e-06 11270
#    19: ll(61.20, -9.26) 1.650e-03 11319  en(10307, 1278686) 6.545e-06 11319
#    20: ll(61.23, -9.27) 1.676e-03 11383  en(10307, 1282514) 6.600e-06 11383
#    21: ll(61.36, -9.30) 1.776e-03 11437  en(10307, 1297037) 6.816e-06 11437
#    22: ll(61.38, -9.30) 1.789e-03 11472  en(10307, 1298889) 6.844e-06 11472
#    23: ll(61.25, -9.39) 1.885e-03 91137  en(4367, 1285831) 7.392e-06 91137
#    24: ll(61.32, -9.40) 1.944e-03 91207  en(4367, 1293568) 7.519e-06 91207
#    25: ll(61.34, -9.41) 1.963e-03 91376  en(4367, 1296061) 7.561e-06 91376
#    26: ll(61.37, -9.41) 1.986e-03 91595  en(4367, 1298908) 7.608e-06 91595
# 131072 total OSGB36
