
# -*- coding: utf-8 -*-

u'''Ordinance Survey Grid References (OSGR) classes L{Osgr} an L{OSGRError}
and functions L{parseOSGR} and L{toOsgr}.

Pure Python implementation of OS Grid Reference functions using an
ellipsoidal earth model, transcribed from JavaScript originals by
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

@newfield example: Example, Examples
'''

from pygeodesy.basics import halfs2, _isnotError, issubclassof, \
                             map1, property_RO
from pygeodesy.datum import Datums
from pygeodesy.dms import parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.fmath import fdot, fpowers, Fsum, fsum_
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import EasNor2Tuple, LatLonDatum3Tuple, \
                           _NamedBase, nameof, _xnamed
from pygeodesy.streprs import enstr2
from pygeodesy.utily import degrees90, degrees180, false2f, sincos2

from math import cos, radians, sin, sqrt, tan

# all public contants, classes and functions
__all__ = _ALL_LAZY.osgr
__version__ = '20.03.23'

_10um    = 1e-5    #: (INTERNAL) 0.01 millimeter (C{meter})
_100km   = 100000  #: (INTERNAL) 100 km (int meter)

_A0, _B0 = radians(49), radians(-2)  #: (INTERNAL) NatGrid true origin, 49°N,2°W.
_E0, _N0 = 400e3, -100e3  #: (INTERNAL) East-/northing of true origin (C{meter}).
_F0      = 0.9996012717   #: (INTERNAL) NatGrid scale of central meridian (C{float}).

_OSGB36  = Datums.OSGB36  #: (INTERNAL) Airy130 ellipsoid
_TRIPS   = 32  #: (INTERNAL) Convergence


def _ll2datum(ll, datum, name):
    '''(INTERNAL) Convert datum if needed.
    '''
    if datum and ll.datum != datum:
        try:
            ll = ll.convertDatum(datum)
        except AttributeError:
            raise TypeError('no %s.convertDatum: %r' % (name, ll))
    return ll


def _M(Mabcd, a):
    '''(INTERNAL) Compute meridional arc.
    '''
    a_ = a - _A0
    _a = a + _A0
    return fdot(Mabcd, a_, -sin(a_)     * cos(_a),
                            sin(a_ * 2) * cos(_a * 2),
                           -sin(a_ * 3) * cos(_a * 3))


class OSGRError(ValueError):
    '''Ordinance Survey Grid References (OSGR) parse or other L{Osgr} issue.
    '''
    pass


class Osgr(_NamedBase):
    '''Ordinance Survey Grid References (OSGR) coordinate.
    '''
    _datum    = _OSGB36  #: (INTERNAL) Datum (L{Datum})
    _easting  = 0        #: (INTERNAL) Easting (C{meter}).
    _latlon   = None     #: (INTERNAL) Cache B{C{_toLatlon}}.
    _northing = 0        #: (INTERNAL) Nothing (C{meter}).

    def __init__(self, easting, northing, name=''):
        '''New L{Osgr} National Grid Reference.

           @arg easting: Easting from OS false easting (C{meter}).
           @arg northing: Northing from from OS false northing (C{meter}).
           @kwarg name: Optional name (C{str}).

           @raise OSGRError: Invalid B{C{easting}} or B{C{northing}}.

           @example:

           >>> from pygeodesy import Osgr
           >>> r = Osgr(651409, 313177)
        '''
        if name:
            self.name = name

        self._easting  = false2f(easting, 'easting',   Error=OSGRError)
        self._northing = false2f(northing, 'northing', Error=OSGRError)

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    def parse(self, strOSGR):
        '''Parse a string to an Osgr instance.

           For more details, see function L{parseOSGR} in this module L{osgr}.
        '''
        return parseOSGR(strOSGR)

    def toLatLon(self, LatLon=None, datum=Datums.WGS84):
        '''Convert this OSGR coordinate to an (ellipsoidal) geodetic
           point.

           I{Note formulation implemented here due to Thomas, Redfearn,
           etc. is as published by OS, but is inferior to Krüger as
           used by e.g. Karney 2011.}

           @kwarg LatLon: Optional ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg datum: Optional datum to use (C{Datum}).

           @return: The geodetic point (B{C{LatLon}}) or a
                    L{LatLonDatum3Tuple}C{(lat, lon, datum)}
                    if B{C{LatLon}} is C{None}.

           @raise OSGRError: No convergence.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal or if
                             B{C{datum}} conversion failed.

           @example:

           >>> from pygeodesy import ellipsoidalVincenty as eV
           >>> g = Osgr(651409.903, 313177.270)
           >>> p = g.toLatLon(eV.LatLon)  # 52°39′28.723″N, 001°42′57.787″E
           >>> # to obtain (historical) OSGB36 lat-/longitude point
           >>> p = g.toLatLon(eV.LatLon, datum=Datums.OSGB36)  # 52°39′27.253″N, 001°43′04.518″E
        '''
        if self._latlon:
            return self._latlon3(LatLon, datum)

        E = self.datum.ellipsoid  # _OSGB36.ellipsoid, Airy130
        a_F0 = E.a * _F0
        b_F0 = E.b * _F0

        e, n = self.easting, self.northing
        n_N0 = n - _N0

        a, m = _A0, n_N0
        sa = Fsum(a)
        for _ in range(_TRIPS):
            a = sa.fsum_(m / a_F0)
            m = n_N0 - b_F0 * _M(E.Mabcd, a)  # meridional arc
            if abs(m) < _10um:
                break
        else:
            t = self.classname, self.toStr(prec=-3), 'toLatLon'
            raise OSGRError('no convergence %s(%s).%s' % t)

        sa, ca = sincos2(a)

        s = E.e2s2(sa)  # r, v = E.roc2_(sa, _F0)
        v = a_F0 / sqrt(s)  # nu
        r = v * E.e12 / s  # rho = a_F0 * E.e12 / pow(s, 1.5) == a_F0 * E.e12 / (s * sqrt(s))

        vr = v / r  # == s / E.e12
        x2 = vr - 1  # η2
        ta = tan(a)

        v3, v5, v7 = fpowers(v, 7, 3)  # PYCHOK false!
        ta2, ta4, ta6 = fpowers(ta**2, 3)  # PYCHOK false!

        tar = ta / r
        V4 = (a,
              tar / (  2 * v),
              tar / ( 24 * v3) * fdot((1, 3, -9), 5 + x2, ta2, ta2 * x2),
              tar / (720 * v5) * fdot((61, 90, 45), 1, ta2, ta4))

        csa = 1.0 / ca
        X5 = (_B0,
              csa / v,
              csa / (   6 * v3) * fsum_(vr, ta, ta),
              csa / ( 120 * v5) * fdot((5, 28, 24), 1, ta2, ta4),
              csa / (5040 * v7) * fdot((61, 662, 1320, 720), 1, ta2, ta4, ta6))

        d, d2, d3, d4, d5, d6, d7 = fpowers(e - _E0, 7)  # PYCHOK false!
        a = fdot(V4, 1,    -d2, d4, -d6)
        b = fdot(X5, 1, d, -d3, d5, -d7)

        self._latlon = _LLEB(degrees90(a), degrees180(b), datum=self.datum, name=self.name)
        return self._latlon3(LatLon, datum)

    def _latlon3(self, LatLon, datum):
        '''(INTERNAL) Convert cached LatLon
        '''
        ll = self._latlon
        if LatLon is None:
            if datum and datum != ll.datum:
                raise TypeError('no %s.convertDatum: %r' % (LatLon, ll))
            return _xnamed(LatLonDatum3Tuple(ll.lat, ll.lon, ll.datum), ll.name)
        elif issubclassof(LatLon, _LLEB):
            ll = _xnamed(LatLon(ll.lat, ll.lon, datum=ll.datum), ll.name)
            return _ll2datum(ll, datum, 'LatLon')
        raise _isnotError(_LLEB.__name__, LatLon=LatLon)

    def toStr(self, prec=10, sep=' '):  # PYCHOK expected
        '''Return a string representation of this OSGR coordinate.

           Note that OSGR coordinates are truncated, not rounded
           (unlike UTM grid references).

           @kwarg prec: Optional number of digits (C{int}).
           @kwarg sep: Optional C{join} separator (C{str}).

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
            return chr(ord('A') + i)

        e, n, s = self._easting, self._northing, ','
        if prec > 0:
            E, e = divmod(e, _100km)
            N, n = divmod(n, _100km)
            E, N = int(E), int(N)
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                return ''
            N = 19 - N
            EN = _i2c( N - (N % 5) + (E + 10) // 5) + \
                 _i2c((N * 5) % 25 + (E % 5))

            t = enstr2(e, n, prec, EN)
            s = sep

        elif -6 < prec < 0:
            w = 6 + 1 - prec
            t = ['%0*.*f' % (w, -prec, t) for t in (e, n)]
        else:
            t = ['%06d' % int(t) for t in (e, n)]
        return s.join(t)

    def toStr2(self, prec=10, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return a string representation of this OSGR coordinate.

           @kwarg prec: Optional number of digits (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator to join (C{str}).

           @return: This OSGR (C{str}) "[G:00B, E:meter, N:meter]" or
                    "OSGR:meter,meter" if B{C{prec}} is non-positive.
        '''
        t = self.toStr(prec=prec, sep=' ')
        if prec > 0:
            t = sep.join('%s:%s' % t for t in zip('GEN', t.split()))
        else:
            t = '%s:%s' % (Osgr.__name__.upper(), t)
        return fmt % (t,)


def parseOSGR(strOSGR, Osgr=Osgr, name=''):
    '''Parse an OSGR coordinate string to an Osgr instance.

       Accepts standard OS Grid References like 'SU 387 148',
       with or without whitespace separators, from 2- up to
       10-digit references (1 m × 1 m square), or fully
       numeric, comma-separated references in metres, for
       example '438700,114800'.

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
        g = ord(G.upper()) - ord('A')
        if g > 7:
            g -= 1
        return g

    def _s2f(g):
        return float(g.strip())

    def _s2i(G, g):
        m = g + '00000'  # std to meter
        return int(str(G) + m[:5])

    s = strOSGR.strip()
    try:
        g = s.split(',')
        if len(g) == 2:  # "easting,northing"
            if len(s) < 13:
                raise ValueError  # caught below
            e, n = map(_s2f, g)

        else:  # "GR easting northing"

            g, s = s[:2], s[2:].strip()

            e, n = map(_c2i, g)
            n, m = divmod(n, 5)
            E = ((e - 2) % 5) * 5 + m
            N = 19 - (e // 5) * 5 - n
            if 0 > E or E > 6 or \
               0 > N or N > 12:
                raise ValueError  # caught below

            g = s.split()
            if len(g) == 1:  # no whitespace
                e, n = halfs2(s)
            elif len(g) == 2:
                e, n = g
            else:
                raise ValueError  # caught below

            e = _s2i(E, e)
            n = _s2i(N, n)

    except ValueError:
        raise OSGRError('%s invalid: %r' % ('strOSGR', strOSGR))

    r = EasNor2Tuple(e, n) if Osgr is None else Osgr(e, n)
    return _xnamed(r, name)


def toOsgr(latlon, lon=None, datum=Datums.WGS84, Osgr=Osgr, name='',
                                               **Osgr_kwds):
    '''Convert a lat-/longitude point to an OSGR coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal) geodetic
                    C{LatLon} point.
       @kwarg lon: Optional longitude in degrees (scalar or C{None}).
       @kwarg datum: Optional datum to convert (C{Datum}).
       @kwarg Osgr: Optional class to return the OSGR coordinate
                    (L{Osgr}) or C{None}.
       @kwarg name: Optional B{C{Osgr}} name (C{str}).
       @kwarg Osgr_kwds: Optional, additional B{C{Osgr}} keyword
                         arguments, ignored if B{C{Osgr=None}}.

       @return: The OSGR coordinate (B{C{Osgr}}) or an
                L{EasNor2Tuple}C{(easting, northing)} if B{C{Osgr}}
                is C{None}.

       @raise TypeError: Non-ellipsoidal B{C{latlon}} or B{C{datum}}
                         conversion failed.

       @raise OSGRError: Invalid B{C{latlon}} or B{C{lon}}.

       @example:

       >>> p = LatLon(52.65798, 1.71605)
       >>> r = toOsgr(p)  # TG 51409 13177
       >>> # for conversion of (historical) OSGB36 lat-/longitude:
       >>> r = toOsgr(52.65757, 1.71791, datum=Datums.OSGB36)
    '''
    if not isinstance(latlon, _LLEB):
        # XXX fix failing _LLEB.convertDatum()
        latlon = _LLEB(*parseDMS2(latlon, lon), datum=datum)
    elif lon is not None:
        raise OSGRError('%s not %s: %r' % ('lon', None, lon))
    elif not name:  # use latlon.name
        name = nameof(latlon)

    E = _OSGB36.ellipsoid

    ll = _ll2datum(latlon, _OSGB36, 'latlon')
    a, b = map1(radians, ll.lat, ll.lon)

    sa, ca = sincos2(a)

    s = E.e2s2(sa)  # r, v = E.roc2_(sa, _F0); r = v / r
    v = E.a * _F0 / sqrt(s)  # nu
    r = s / E.e12  # nu / rho == v / (v * E.e12 / s) == s / E.e12

    x2 = r - 1  # η2
    ta = tan(a)

    ca3, ca5 = fpowers(ca, 5, 3)  # PYCHOK false!
    ta2, ta4 = fpowers(ta, 4, 2)  # PYCHOK false!

    vsa = v * sa
    I4 = (E.b * _F0 * _M(E.Mabcd, a) + _N0,
         (vsa /   2) * ca,
         (vsa /  24) * ca3 * fsum_(5, -ta2, 9 * x2),
         (vsa / 720) * ca5 * fsum_(61, ta4, -58 * ta2))

    V4 = (_E0,
         (v        * ca),
         (v /   6) * ca3 * (r - ta2),
         (v / 120) * ca5 * fdot((-18, 1, 14, -58), ta2, 5 + ta4, x2, ta2 * x2))

    d, d2, d3, d4, d5, d6 = fpowers(b - _B0, 6)  # PYCHOK false!
    n = fdot(I4, 1, d2, d4, d6)
    e = fdot(V4, 1, d,  d3, d5)

    if Osgr is None:
        r = EasNor2Tuple(e, n)
    else:
        r = Osgr(e, n, **Osgr_kwds)
        if lon is None and isinstance(latlon, _LLEB):
            r._latlon = latlon  # XXX weakref(latlon)?
    return _xnamed(r, name)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
