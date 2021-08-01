
# -*- coding: utf-8 -*-

u'''Ordinance Survey Grid References (OSGR) references.

Classes L{Osgr} and L{OSGRError} and functions L{parseOSGR} and L{toOsgr}.

Pure Python implementation of OS Grid Reference functions using an
ellipsoidal earth model, transcoded from JavaScript originals by
I{(C) Chris Veness 2005-2016} published under the same MIT Licence**, see
U{OS National Grid<https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>}
and U{Module osgridref
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-osgridref.html>}.

OSGR provides geocoordinate references for UK mapping purposes, converted
in 2015 to work with WGS84 datum by default or OSGB36 as option.

See U{Guide<https://www.OrdnanceSurvey.co.UK/docs/support/guide-coordinate-systems-great-britain.pdf>},
U{Proposed Changes<https://www.OrdnanceSurvey.co.UK/blog/2014/09/proposed-changes-to-latitude-and-longitude-representation-on-paper-maps-tell-us-your-thoughts>},
U{Confirmation<https://www.OrdnanceSurvey.co.UK/blog/2014/12/confirmation-on-changes-to-latitude-and-longitude>}
and U{Ordnance Survey National Grid<https://WikiPedia.org/wiki/Ordnance_Survey_National_Grid>}.

See also Karney U{'Transverse Mercator with an accuracy of a few nanometers'
<https://Arxiv.org/pdf/1002.1417v3.pdf>}, 2011 (building on Krüger
U{'Konforme Abbildung des Erdellipsoids in der Ebene'
<https://bib.GFZ-Potsdam.DE/pub/digi/krueger2.pdf>}, 1912), Seidel
U{'Die Mathematik der Gauß-Krueger-Abbildung'
<https://Henrik-Seidel.GMXhome.DE/gausskrueger.pdf>}, 2006 and
U{Transverse Mercator: Redfearn series
<https://WikiPedia.org/wiki/Transverse_Mercator:_Redfearn_series>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import halfs2, map1, _xsubclassof
from pygeodesy.datums import Datums, _ellipsoidal_datum, _WGS84
from pygeodesy.dms import parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _parseX, _TypeError, _ValueError
from pygeodesy.fmath import fdot, fpowers, Fsum, fsum_
from pygeodesy.interns import NN, _A_, _COMMA_, _COMMASPACE_, _DOT_, \
                             _convergence_, _float, _latlon_, _no_, \
                             _not_, _SPACE_, _1_0, _2_0, _6_0, \
                             _24_0, _120_0, _720_0
from pygeodesy.interns import _COLON_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedBase, nameof, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, LatLonDatum3Tuple
from pygeodesy.props import Property_RO, property_RO
from pygeodesy.streprs import enstr2, Fmt, _xzipairs, _0wd, _0wpF
from pygeodesy.units import Easting, Lam_, Northing, Phi_, Scalar, \
                           _10um, _100km
from pygeodesy.utily import degrees90, degrees180, sincos2

from math import cos, radians, sin, sqrt, tan

__all__ = _ALL_LAZY.osgr
__version__ = '21.07.31'

_100_000 =  int(_100km)  # 100 km (int C{meter})
_5040_0  = _float(5040)

_A0 = Phi_(49)              # NatGrid true origin latitude, 49°N
_B0 = Lam_(-2)              # NatGrid true origin longitude, 2°W
_E0 = Easting(400e3)        # Easting of true origin (C{meter})
_N0 = Northing(-_100km)     # Northing of true origin (C{meter})
_F0 = Scalar(0.9996012717)  # NatGrid scale of central meridian (C{float})

_OSGB36      =  Datums.OSGB36  # Airy130 ellipsoid
_no_toDatum_ = 'no .toDatum'
_ord_A       =  ord(_A_)
_TRIPS       =  33  # .toLatLon convergence


def _ll2datum(ll, datum, name):
    '''(INTERNAL) Convert datum if needed.
    '''
    if datum not in (None, ll.datum):
        try:
            ll = ll.toDatum(datum)
        except AttributeError:
            raise _TypeError(name, ll, txt=Fmt.PAREN(_no_toDatum_, datum.name))
    return ll


def _M(Mabcd, a):
    '''(INTERNAL) Compute meridional arc.
    '''
    a_ = a - _A0
    _a = a + _A0
    return fdot(Mabcd, a_, -sin(a_)     * cos(_a),
                            sin(a_ * 2) * cos(_a * 2),
                           -sin(a_ * 3) * cos(_a * 3))


class OSGRError(_ValueError):
    '''Ordinance Survey Grid References (OSGR) parse or other L{Osgr} issue.
    '''
    pass


class Osgr(_NamedBase):
    '''Ordinance Survey Grid References (OSGR) coordinate.
    '''
    _datum     = _OSGB36  # default datum (L{Datum})
    _easting   =  0       # Easting (C{meter})
    _iteration =  None    # iteration number (C{int})
    _latlon    =  None    # cached B{C{_toLatlon}}
    _northing  =  0       # Nothing (C{meter})

    def __init__(self, easting, northing, datum=None, name=NN):
        '''New L{Osgr} National Grid Reference.

           @arg easting: Easting from OS false easting (C{meter}).
           @arg northing: Northing from from OS false northing (C{meter}).
           @kwarg datum: Default datum (C{Datums.OSGB36}).
           @kwarg name: Optional name (C{str}).

           @raise OSGRError: Invalid or negative B{C{easting}} or
                             B{C{northing}} or B{C{datum}} not
                             C{Datums.OSBG36}.

           @example:

            >>> from pygeodesy import Osgr
            >>> r = Osgr(651409, 313177)
        '''
        self._easting  = Easting( easting,  Error=OSGRError, osgr=True)
        self._northing = Northing(northing, Error=OSGRError, osgr=True)

        if datum not in (None, _OSGB36):
            try:
                if _ellipsoidal_datum(datum) != _OSGB36:
                    raise ValueError
            except (TypeError, ValueError):
                raise OSGRError(datum=datum)
        if name:
            self.name = name

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

    @property_RO
    def iteration(self):
        '''Get the most recent C{Osgr.toLatLon} iteration number
           (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    def parse(self, strOSGR, name=NN):
        '''Parse a string to a similar L{Osgr} instance.

           @arg strOSGR: The OSGR reference (C{str}),
                         see function L{parseOSGR}.
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar instance (L{Osgr})

           @raise OSGRError: Invalid B{C{strOSGR}}.
        '''
        return parseOSGR(strOSGR, Osgr=self.classof, name=name or self.name)

    def toLatLon(self, LatLon=None, datum=_WGS84):
        '''Convert this OSGR coordinate to an (ellipsoidal) geodetic
           point.

           While OS grid references are based on the OSGB36 datum, the
           I{Ordnance Survey} have deprecated the use of OSGB36 for
           lat-/longitude coordinates (in favour of WGS84). Hence, this
           method returns WGS84 by default with OSGB36 as an option,
           U{see<https://www.OrdnanceSurvey.co.UK/blog/2014/12/2>}.

           I{Note formulation implemented here due to Thomas, Redfearn,
           etc. is as published by OS, but is inferior to Krüger as
           used by e.g. Karney 2011.}

           @kwarg LatLon: Optional ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg datum: Optional datum to convert to (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{Ellipsoid2}
                         or L{a_f2Tuple}).

           @return: The geodetic point (B{C{LatLon}}) or a
                    L{LatLonDatum3Tuple}C{(lat, lon, datum)}
                    if B{C{LatLon}} is C{None}.

           @raise OSGRError: No convergence.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal or
                             B{C{datum}} is invalid or conversion failed.

           @example:

            >>> from pygeodesy import ellipsoidalVincenty as eV
            >>> g = Osgr(651409.903, 313177.270)
            >>> p = g.toLatLon(eV.LatLon)  # 52°39′28.723″N, 001°42′57.787″E
            >>> # to obtain (historical) OSGB36 lat-/longitude point
            >>> p = g.toLatLon(eV.LatLon, datum=Datums.OSGB36)  # 52°39′27.253″N, 001°43′04.518″E
        '''
        if self._latlon:
            return self._latlon3(LatLon, datum)

        E = self.datum.ellipsoid  # _Datums_OSGB36.ellipsoid, Airy130
        a_F0 = E.a * _F0
        b_F0 = E.b * _F0

        e, n = self.easting, self.northing
        n_N0 = n - _N0

        a, m = _A0, n_N0
        A = Fsum(a)
        for self._iteration in range(1, _TRIPS):
            a = A.fsum_(m / a_F0)
            m = n_N0 - b_F0 * _M(E.Mabcd, a)  # meridional arc
            if abs(m) < _10um:
                break
        else:
            t = _DOT_(Fmt.PAREN(self.classname, self.toStr(prec=-3)),
                                self.toLatLon.__name__)
            raise OSGRError(_no_(_convergence_), txt=t)
        sa, ca = sincos2(a)

        s = E.e2s2(sa)  # r, v = E.roc2_(sa, _F0)
        v = a_F0 / sqrt(s)  # nu
        r = v * E.e12 / s  # rho = a_F0 * E.e12 / pow(s, 1.5) == a_F0 * E.e12 / (s * sqrt(s))

        vr = v / r  # == s / E.e12
        x2 = vr - _1_0  # η2
        ta = tan(a)

        v3, v5, v7 = fpowers(v, 7, alts=3)  # PYCHOK false!
        ta2, ta4, ta6 = fpowers(ta**2, 3)  # PYCHOK false!

        d1, d2, d3, d4, d5, d6, d7 = fpowers(e - _E0, 7)  # PYCHOK false!

        t = ta / r
        a = fsum_(a,
                 -d2 * t / (  _2_0 * v),
                  d4 * t / ( _24_0 * v3) * fsum_(5, x2, 3 * ta2, -9 * ta2 * x2),
                 -d6 * t / (_720_0 * v5) * fsum_(61, 90 * ta2, 45 * ta4))

        t = _1_0 / ca
        b = fsum_(_B0,
                   d1 * t /            v,
                  -d3 * t / (   _6_0 * v3) * fsum_(vr, ta2, ta2),
                   d5 * t / ( _120_0 * v5) * fsum_(5,   28 * ta2,   24 * ta4),
                  -d7 * t / (_5040_0 * v7) * fsum_(61, 662 * ta2, 1320 * ta4, 720 * ta6))

        r = _LLEB(degrees90(a), degrees180(b), datum=self.datum, name=self.name)
        r._iteration = self._iteration  # only ellipsoidal LatLon
        self._latlon = r
        return self._latlon3(LatLon, datum)

    def _latlon3(self, LatLon, datum):
        '''(INTERNAL) Convert cached latlon to C{LatLon}
        '''
        ll = self._latlon
        if LatLon is None:
            r = _ll2datum(ll, datum, LatLonDatum3Tuple.__name__)
            r =  LatLonDatum3Tuple(r.lat, r.lon, r.datum)
        else:  # must be ellipsoidal
            _xsubclassof(_LLEB, LatLon=LatLon)
            r = _ll2datum(ll, datum, LatLon.__name__)
            r =  LatLon(r.lat, r.lon, datum=r.datum)
        r._iteration = ll._iteration
        return _xnamed(r, nameof(ll))

    def toRepr(self, prec=10, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this OSGR coordinate.

           @kwarg prec: Optional number of digits (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator to join (C{str}).

           @return: This OSGR (C{str}) "[G:00B, E:meter, N:meter]" or
                    "[OSGR:meter,meter]" if B{C{prec}} is non-positive.
        '''
        t = self.toStr(prec=prec, sep=None)
        if prec > 0:
            t = _xzipairs('GEN', t, sep=sep, fmt=fmt)
        else:
            t = _COLON_(Osgr.__name__.upper(), t)
            if fmt:
                t = fmt % (t,)
        return t

    def toStr(self, prec=10, sep=_SPACE_):  # PYCHOK expected
        '''Return a string representation of this OSGR coordinate.

           Note that OSGR coordinates are truncated, not rounded
           (unlike UTM grid references).

           @kwarg prec: Optional number of digits (C{int}).
           @kwarg sep: Optional C{join} separator (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.

           @return: This OSGR as C{"EN easting northing"} or as
                    C{"easting,northing"} if B{C{prec}} is non-positive
                    (C{str}).

           @raise ValueError: Invalid B{C{prec}}.

           @example:

            >>> r = Osgr(651409, 313177)
            >>> str(r)  # TG 5140 1317
            >>> r.toStr(prec=0)  # 651409,313177
        '''
        def _i2c(i):
            if i > 7:
                i += 1
            return chr(_ord_A + i)

        e, n, s = self._easting, self._northing, _COMMA_
        if prec > 0:
            E, e = divmod(e, _100_000)
            N, n = divmod(n, _100_000)
            E, N = int(E), int(N)
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                return NN
            N = 19 - N
            EN = _i2c( N - (N % 5) + (E + 10) // 5) + \
                 _i2c((N * 5) % 25 + (E % 5))

            t = enstr2(e, n, prec, EN)
            s = sep

        elif -6 < prec < 0:
            w = 6 + 1 - prec
            t = [_0wpF(w, -prec, t) for t in (e, n)]
        else:
            t = [_0wd(6, int(t)) for t in (e, n)]

        return tuple(t) if s is None else s.join(t)


def parseOSGR(strOSGR, Osgr=Osgr, name=NN):
    '''Parse a string representing an OSGR grid reference,
       consisting of C{"[grid] easting northing"}.

       Accepts standard OS Grid References like 'SU 387 148',
       with or without whitespace separators, from 2- up to
       10-digit references (1 m × 1 m square), or all-numeric,
       comma-separated references in meters, for example
       '438700,114800'.

       @arg strOSGR: An OSGR coordinate (C{str}).
       @kwarg Osgr: Optional class to return the OSGR
                    coordinate (L{Osgr}) or C{None}.
       @kwarg name: Optional B{C{Osgr}} name (C{str}).

       @return: The OSGR coordinate (B{C{Osgr}}) or an
                L{EasNor2Tuple}C{(easting, northing)} if
                B{C{Osgr}} is C{None}.

       @raise OSGRError: Invalid B{C{strOSGR}}.

       @example:

        >>> g = parseOSGR('TG 51409 13177')
        >>> str(g)  # TG 51409 13177
        >>> g = parseOSGR('TG5140913177')
        >>> str(g)  # TG 51409 13177
        >>> g = parseOSGR('TG51409 13177')
        >>> str(g)  # TG 51409 13177
        >>> g = parseOSGR('651409,313177')
        >>> str(g)  # TG 51409 13177
        >>> g.toStr(prec=0)  # 651409,313177
    '''
    def _c2i(G):
        g = ord(G.upper()) - _ord_A
        if g > 7:
            g -= 1
        return g

    def _s2f(g):
        return float(g.strip())

    def _s2i(G, g):
        m = g + '00000'  # std to meter
        return int(str(G) + m[:5])

    def _OSGR_(strOSGR, Osgr, name):
        s = strOSGR.strip()
        g = s.split(_COMMA_)
        if len(g) == 2:  # "easting,northing"
            if len(s) < 13:
                raise ValueError
            e, n = map(_s2f, g)

        else:  # "GR easting northing"
            g, s = s[:2], s[2:].strip()

            e, n = map(_c2i, g)
            n, m = divmod(n, 5)
            E = ((e - 2) % 5) * 5 + m
            N = 19 - (e // 5) * 5 - n
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                raise ValueError

            g = s.split()
            if len(g) == 1:  # no whitespace
                e, n = halfs2(s)
            elif len(g) == 2:
                e, n = g
            else:
                raise ValueError

            e = _s2i(E, e)
            n = _s2i(N, n)

        return EasNor2Tuple(e, n, name=name) if Osgr is None else \
               _xnamed(Osgr(e, n), name)

    return _parseX(_OSGR_, strOSGR, Osgr, name,
                           strOSGR=strOSGR, Error=OSGRError)


def toOsgr(latlon, lon=None, datum=_WGS84, Osgr=Osgr, name=NN,
                                         **Osgr_kwds):
    '''Convert a lat-/longitude point to an OSGR coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal) geodetic
                    C{LatLon} point.
       @kwarg lon: Optional longitude in degrees (scalar or C{None}).
       @kwarg datum: Optional datum to convert B{C{lat, lon}} from
                     (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                     L{a_f2Tuple}).
       @kwarg Osgr: Optional class to return the OSGR coordinate
                    (L{Osgr}) or C{None}.
       @kwarg name: Optional B{C{Osgr}} name (C{str}).
       @kwarg Osgr_kwds: Optional, additional B{C{Osgr}} keyword
                         arguments, ignored if C{B{Osgr} is None}.

       @return: The OSGR coordinate (B{C{Osgr}}) or an
                L{EasNor2Tuple}C{(easting, northing)} if B{C{Osgr}}
                is C{None}.

       @raise OSGRError: Invalid B{C{latlon}} or B{C{lon}}.

       @raise TypeError: Non-ellipsoidal B{C{latlon}} or invalid
                         B{C{datum}} or conversion failed.

       @example:

        >>> p = LatLon(52.65798, 1.71605)
        >>> r = toOsgr(p)  # TG 51409 13177
        >>> # for conversion of (historical) OSGB36 lat-/longitude:
        >>> r = toOsgr(52.65757, 1.71791, datum=Datums.OSGB36)
    '''
    if not isinstance(latlon, _LLEB):
        # XXX fix failing _LLEB.toDatum()
        latlon = _LLEB(*parseDMS2(latlon, lon), datum=datum)
    elif lon is not None:
        raise OSGRError(lon=lon, txt=_not_(None))
    elif not name:  # use latlon.name
        name = nameof(latlon)

    # if necessary, convert to OSGB36 first
    ll = _ll2datum(latlon, _OSGB36, _latlon_)
    try:
        a, b = ll.philam
    except AttributeError:
        a, b = map1(radians, ll.lat, ll.lon)
    sa, ca = sincos2(a)

    E = _OSGB36.ellipsoid

    s = E.e2s2(sa)  # r, v = E.roc2_(sa, _F0); r = v / r
    v = E.a * _F0 / sqrt(s)  # nu
    r = s / E.e12  # nu / rho == v / (v * E.e12 / s) == s / E.e12

    x2 = r - _1_0  # η2
    ta = tan(a)

    ca3, ca5 = fpowers(ca, 5, alts=3)  # PYCHOK false!
    ta2, ta4 = fpowers(ta, 4, alts=2)  # PYCHOK false!

    d1, d2, d3, d4, d5, d6 = fpowers(b - _B0, 6)  # PYCHOK false!

    t = fsum_(-18 * ta2, 5, ta4, 14 * x2, -58 * ta2 * x2)
    e = fsum_(_E0,
               d1 * v          * ca,
               d3 * v /   _6_0 * ca3 * (r - ta2),
               d5 * v / _120_0 * ca5 * t)

    t = v * sa
    n = fsum_(_N0,
              _F0 * E.b * _M(E.Mabcd, a),
               d2 * t /   _2_0 * ca,
               d4 * t /  _24_0 * ca3 * fsum_(5, -ta2,   9 * x2),
               d6 * t / _720_0 * ca5 * fsum_(61, ta4, -58 * ta2))

    if Osgr is None:
        r = EasNor2Tuple(e, n)
    else:
        r = Osgr(e, n, datum=_OSGB36, **Osgr_kwds)
        if lon is None and isinstance(latlon, _LLEB):
            r._latlon = latlon  # XXX weakref(latlon)?
    return _xnamed(r, name or nameof(latlon))

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
