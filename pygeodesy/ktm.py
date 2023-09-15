
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ class U{TransverseMercator
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}
based on I{Krüger} series.  See also I{Karney}'s utility U{TransverseMercatorProj
<https://GeographicLib.SourceForge.io/C++/doc/TransverseMercatorProj.1.html>}.

Following and further below is a copy of I{Karney}'s U{TransverseMercator.hpp
<https://GeographicLib.SourceForge.io/C++/doc/TransverseMercator_8hpp_source.html>}
file C{Header}.

This implementation follows closely JHS 154, ETRS89 - I{järjestelmään liittyvät
karttaprojektiot, tasokoordinaatistot ja karttalehtijako} (Map projections, plane
coordinates, and map sheet index for ETRS89), published by JUHTA, Finnish Geodetic
Institute, and the National Land Survey of Finland (2006).  The relevant section
is available as the U{2008 PDF file
<http://Docs.JHS-suositukset.FI/jhs-suositukset/JHS154/JHS154_liite1.pdf>}.

This is a straight transcription of the formulas in this paper with the
following exceptions:

 - Use of 6th order series instead of 4th order series. This reduces the
   error to about 5 nm for the UTM range of coordinates (instead of 200 nm),
   with a speed penalty of only 1%,

 - Use Newton's method instead of plain iteration to solve for latitude
   in terms of isometric latitude in the Reverse method,

 - Use of Horner's representation for evaluating polynomials and Clenshaw's
   method for summing trigonometric series,

 - Several modifications of the formulas to improve the numerical accuracy,

 - Evaluating the convergence and scale using the expression for the
   projection or its inverse.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2023)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, isodd, neg, neg_, _reverange
from pygeodesy.constants import INF, _K0_UTM, PI, PI_2, _0_0s, _0_0, \
                               _1_0, _90_0, _copysignINF
# from pygeodesy.datums import _spherical_datum  # _MODS
# from pygeodesy.ellipsoids import _EWGS84  # from .karney
from pygeodesy.errors import _ValueError, _xkwds_get, _Xorder
from pygeodesy.fmath import hypot, hypot1
from pygeodesy.fsums import fsum1f_
from pygeodesy.interns import NN, _COMMASPACE_, _singular_
from pygeodesy.karney import _atan2d, _diff182, _fix90, _polynomial, \
                             _norm180, _unsigned2,  _EWGS84, _NamedBase
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS,  _pairs
# from pygeodesy.named import _NamedBase  # from .karney
from pygeodesy.namedTuples import Forward4Tuple, Reverse4Tuple
from pygeodesy.props import property_doc_, Property, Property_RO, \
                           _update_all
# from pygeodesy.streprs import pairs as _pairs  # from .lazily
from pygeodesy.units import Degrees, Scalar_, _1mm as _TOL_10  # PYCHOK used!
from pygeodesy.utily import atand, _loneg, sincos2, sincos2d_

from cmath import polar
from math import atan2, asinh, cos, cosh, degrees, fabs, sin, sinh, sqrt, tanh

__all__ = _ALL_LAZY.ktm
__version__ = '23.09.07'


class KTMError(_ValueError):
    '''Error raised for L{KTransverseMercator} and L{KTransverseMercator.forward} issues.
    '''
    pass


class KTransverseMercator(_NamedBase):
    '''I{Karney}'s C++ class U{TransverseMercator<https://GeographicLib.SourceForge.io/
       C++/doc/classGeographicLib_1_1TransverseMercator.html>} transcoded to pure
       Python, following is a partial copy of I{Karney}'s documentation.

       Transverse Mercator projection based on Krüger's method which evaluates the
       projection and its inverse in terms of a series.

       There's a singularity in the projection at I{phi = 0, lam - lam0 = +/- (1 - e)
       90}, about +/- 82.6 degrees for WGS84, where I{e} is the eccentricity.  Beyond
       this point, the series ceases to converge and the results from this method
       will be garbage.  I{To be on the safe side, don't use this method if the
       angular distance from the central meridian exceeds (1 - 2e) x 90}, about 75
       degrees for the WGS84 ellipsoid.

       Class L{ExactTransverseMercator} is an alternative implementation of the
       projection using I{exact} formulas which yield accurate (to 8 nm) results
       over the entire ellipsoid.

       The ellipsoid parameters and the central scale are set in the constructor.
       The central meridian (which is a trivial shift of the longitude) is specified
       as the C{lon0} keyword argument of the L{KTransverseMercator.forward} and
       L{KTransverseMercator.reverse} methods.  The latitude of origin is taken to
       be the equator.  There is no provision in this class for specifying a false
       easting or false northing or a different latitude of origin.  However these
       are can be simply included by the calling function.

       The L{KTransverseMercator.forward} and L{KTransverseMercator.reverse} methods
       also return the meridian convergence C{gamma} and scale C{k}.  The meridian
       convergence is the bearing of grid North, the C{y axis}, measured clockwise
       from true North.
    '''
    _E      = _EWGS84
    _k0     = _K0_UTM  # central scale factor
    _lon0   = _0_0     # central meridian
    _mTM    =  6
    _raiser =  False   # throw Error

    def __init__(self, a_earth=_EWGS84, f=None, lon0=0, k0=_K0_UTM, name=NN,
                                                raiser=False, **TMorder):
        '''New L{KTransverseMercator}.

           @kwarg a_earth: This rhumb's earth (L{Ellipsoid}, L{Ellipsoid2},
                           L{a_f2Tuple}, L{Datum}, 2-tuple (C{a, f})) or the
                           equatorial radius (C{scalar}, C{meter}).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     a C{scalar}, ignored otherwise.
           @kwarg lon0: The central meridian (C{degrees180}).
           @kwarg k0: Central scale factor (C{scalar}).
           @kwarg name: Optional name (C{str}).
           @kwarg raiser: If C{True}, throw a L{KTMError} for C{forward}
                          singularities (C{bool}).
           @kwarg TMorder: Keyword argument B{C{TMorder}}, see property C{TMorder}.

           @raise KTMError: Invalid B{C{a_earth}}, B{C{f}} or B{C{TMorder}}.
        '''
        if f is not None:
            self.ellipsoid = a_earth, f
        elif a_earth not in (_EWGS84, None):
            self.ellipsoid = a_earth
        self.lon0 = lon0
        self.k0 = k0
        if name:  # PYCHOK no cover
            self.name = name
        if raiser:
            self.raiser = True
        if TMorder:
            self.TMorder = _xkwds_get(TMorder, TMorder=self._mTM)

    @Property_RO
    def _Alp(self):
        return _Xs(_AlpCoeffs, self.TMorder, self.ellipsoid)

    @Property_RO
    def _b1(self):
        n = self.ellipsoid.n
        if n:  # isEllipsoidal
            m  =  self.TMorder // 2
            B1 = _B1Coeffs[m]
            m +=  1
            b1 = _polynomial(n**2, B1, 0, m) / (B1[m] * (n + _1_0))
        else:  # isSpherical
            b1 = _1_0  # B1[m - 1] / B1[m1] == 1, always
        return b1

    @Property_RO
    def _Bet(self):
        C = _Xs(_BetCoeffs, self.TMorder, self.ellipsoid)
        return tuple(map(neg, C)) if self.f else C  # negated if isEllipsoidal

    @Property
    def ellipsoid(self):
        '''Get the ellipsoid (L{Ellipsoid}).
        '''
        return self._E

    @ellipsoid.setter  # PYCHOK setter!
    def ellipsoid(self, a_earth_f):
        '''Set this rhumb's ellipsoid (L{Ellipsoid}, L{Ellipsoid2}, L{Datum},
           L{a_f2Tuple} or 2-tuple C{(a, f)}).
        '''
        E = _MODS.datums._spherical_datum(a_earth_f, Error=KTMError).ellipsoid
        if self._E != E:
            _update_all(self)
            self._E = E

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

    def forward(self, lat, lon, lon0=None, name=NN):
        '''Forward projection, from geographic to transverse Mercator.

           @arg lat: Latitude of point (C{degrees90}).
           @arg lon: Longitude of point (C{degrees180}).
           @arg lon0: Central meridian of the projection (C{degrees180}).
           @kwarg name: Optional name (C{str}).

           @return: L{Forward4Tuple}C{(easting, northing, gamma, scale)}
                    with C{easting} and C{northing} in C{meter}, unfalsed, the
                    meridian convergence C{gamma} at point in C{degrees180}
                    and the C{scale} of projection at point C{scalar}.  Any
                    value may be C{NAN}, C{NINF} or C{INF} for singularities.

           @raise KTMError: For singularities, iff property C{raiser} is
                            C{True}.
        '''
        lat, _lat = _unsigned2(_fix90(lat))
        lon, _    = _diff182((self.lon0 if lon0 is None else lon0), lon)
        lon, _lon = _unsigned2(lon)
        backside  =  lon > 90
        if backside:  # PYCHOK no cover
            lon = _loneg(lon)
            if lat == 0:
                _lat = True

        sphi, cphi, slam, clam = sincos2d_(lat, lon)
        E = self.ellipsoid
        if cphi and lat != 90:
            t  = sphi / cphi
            tp = E.es_taupf(t)
            h  = hypot(tp, clam)
            if h:
                xip  = atan2(tp, clam)
                etap = asinh(slam / h)  # atanh(sin(lam) / cosh(psi))
                g = _atan2d(slam * tp, clam * hypot1(tp))  # Krueger p 22 (44)
                k =  sqrt(cphi**2 * E.e2 + E.e21) * hypot1(t) / h
            elif self.raiser:
                raise KTMError(lat=lat, lon=lon, lon0=lon0, txt=_singular_)
            else:  # PYCHOK no cover
                xip, etap = _0_0, _copysignINF(slam)
                g, k = copysign0(_90_0, slam), INF
        else:  # PYCHOK no cover
            xip, etap = PI_2, _0_0
            g, k = lon, E.es_c
        y, x, d, t = _Cyxgk4(E, xip, etap, self._Alp)
        g -= d
        k *= t * self._k0_b1

        if backside:  # PYCHOK no cover
            y, g = (PI - y), _loneg(g)
        y *= self._k0_a1
        x *= self._k0_a1
        if _lat:
            y, g = neg_(y, g)
        if _lon:
            x, g = neg_(x, g)

        return Forward4Tuple(x, y, _norm180(g), k, name=name or self.name)

    @property_doc_(''' the central scale factor (C{float}).''')
    def k0(self):
        '''Get the central scale factor (C{float}), aka I{C{scale0}}.
        '''
        return self._k0  # aka scale0

    @k0.setter  # PYCHOK setter!
    def k0(self, k0):
        '''Set the central scale factor (C{float}), aka I{C{scale0}}.

           @raise KTMError: Invalid B{C{k0}}.
        '''
        k0 = Scalar_(k0=k0, Error=KTMError, low=_TOL_10, high=_1_0)
        if self._k0 != k0:  # PYCHOK no cover
            KTransverseMercator._k0_a1._update(self)  # redo ._k0_a1
            KTransverseMercator._k0_b1._update(self)  # redo ._k0_b1
            self._k0 = k0

    @Property_RO
    def _k0_a1(self):
        '''(INTERNAL) Cache C{k0 * _b1 * equatoradius}.
        '''
        return self._k0_b1 * self.equatoradius

    @Property_RO
    def _k0_b1(self):
        '''(INTERNAL) Cache C{k0 * _b1}.
        '''
        return self.k0 * self._b1

    @property_doc_(''' the central meridian (C{degrees180}).''')
    def lon0(self):
        '''Get the central meridian (C{degrees180}).
        '''
        return self._lon0

    @lon0.setter  # PYCHOK setter!
    def lon0(self, lon0):
        '''Set the central meridian (C{degrees180}).

           @raise KTMError: Invalid B{C{lon0}}.
        '''
        self._lon0 = _norm180(Degrees(lon0=lon0, Error=KTMError))

    @property_doc_(''' raise a L{KTMError} for C{forward} singularities (C{bool}).''')
    def raiser(self):
        '''Get the error setting (C{bool}).
        '''
        return self._raiser

    @raiser.setter  # PYCHOK setter!
    def raiser(self, raiser):
        '''Set the error setting (C{bool}), to C{True} to throw a L{KTMError}
           for C{forward} singularities.
        '''
        self._raiser = bool(raiser)

    def reverse(self, x, y, lon0=None, name=NN):
        '''Reverse projection, from transverse Mercator to geographic.

           @arg x: Easting of point (C{meter}).
           @arg y: Northing of point (C{meter}).
           @arg lon0: Central meridian of the projection (C{degrees180}).

           @return: L{Reverse4Tuple}C{(lat, lon, gamma, scale)} with
                    C{lat}- and C{lon}gitude in C{degrees}, I{unfalsed}.
        '''
        eta, _lon = _unsigned2(x / self._k0_a1)
        xi,  _lat = _unsigned2(y / self._k0_a1)
        backside  =  xi > PI_2
        if backside:  # PYCHOK no cover
            xi = PI - xi

        E = self.ellipsoid
        xip, etap, g, k = _Cyxgk4(E, xi, eta, self._Bet)
        t = self._k0_b1
        k = (t / k) if k else _copysignINF(t)  # _over(t, k)
        h, c = sinh(etap), cos(xip)
        if c > 0:
            r =  hypot(h, c)
        else:  # PYCHOK no cover
            r =  fabs(h)
            c = _0_0
        if r:
            lon = _atan2d(h, c)  # Krueger p 17 (25)
            s   =  sin(xip)  # Newton for tau
            t   =  E.es_tauf(s / r)
            lat =  atand(t)
            g  += _atan2d(s * tanh(etap), c)  # Krueger p 19 (31)
            k  *=  sqrt(E.e2 / (t**2 + _1_0) + E.e21) * hypot1(t) * r
        else:  # PYCHOK no cover
            lat, lon = _90_0, _0_0
            k *= E.es_c

        if backside:  # PYCHOK no cover
            lon, g = _loneg(lon), _loneg(g)
        if _lat:
            lat, g = neg_(lat, g)
        if _lon:
            lon, g = neg_(lon, g)

        lon += self.lon0 if lon0 is None else _norm180(lon0)
        return Reverse4Tuple(lat, _norm180(lon), _norm180(g), k,
                                   name=name or self.name)

    @Property
    def TMorder(self):
        '''Get the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mTM

    @TMorder.setter  # PYCHOK setter!
    def TMorder(self, order):
        '''Set the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        m = _Xorder(_AlpCoeffs, KTMError, TMorder=order)
        if self._mTM != m:
            _update_all(self)
            self._mTM = m

    def toStr(self, **kwds):
        '''Return a C{str} representation.

           @arg kwds: Optional, overriding keyword arguments.
        '''
        d = dict(ellipsoid=self.ellipsoid, k0=self.k0, TMorder=self.TMorder)
        if self.name:  # PYCHOK no cover
            d.update(name=self.name)
        return _COMMASPACE_.join(_pairs(d, **kwds))


def _cma(a, b0, b1, Cn):
    '''(INTERNAL) Compute complex M{a * b0 - b1 + Cn} with complex
       C{a}, C{b0} and C{b1} and scalar C{Cn}.

       @see: CPython function U{_Py_c_prod<https://GitHub.com/python/
             cpython/blob/main/Objects/complexobject.c>}.

       @note: Python function C{cmath.fsum} is no longer available,
              but stil mentioned in Note 4 of the comments before
              CPython function U{math_fsum<https://GitHub.com/python/
              cpython/blob/main/Modules/mathmodule.c>}
    '''
    r = fsum1f_(a.real * b0.real, -a.imag * b0.imag, -b1.real, Cn)
    j = fsum1f_(a.real * b0.imag,  a.imag * b0.real, -b1.imag)
    return complex(r, j)


def _Cyxgk4(E, xi_, eta_, C):
    '''(INTERNAL) Complex Clenshaw summation with C{B{C}=._Alp}
       or C{B{C}=-._Bet}.
    '''
    x = complex(xi_, eta_)
    if E.f:  # isEllipsoidal
        s,  c  =  sincos2(  xi_  * 2)
        sh, ch = _sinhcosh2(eta_ * 2)
        n = -s
        s = complex(s * ch, c * sh)  # sin(zeta * 2)
        c = complex(c * ch, n * sh)  # cos(zeta * 2)
        a = c * 2  # cos(zeta * 2) * 2

        y0 = y1 = \
        z0 = z1 = complex(0)  # 0+0j
        n  = len(C) - 1  # == .TMorder
        if isodd(n):
            Cn = C[n]
            y0 = complex(Cn)  # +0j
            z0 = complex(Cn * (n * 2))
            n -= 1
        _c = _cma
        while n > 0:
            Cn =  C[n]
            y1 = _c(a, y0, y1, Cn)
            z1 = _c(a, z0, z1, Cn * (n * 2))
            n -=  1
            Cn =  C[n]
            y0 = _c(a, y1, y0, Cn)
            z0 = _c(a, z1, z0, Cn * (n * 2))
            n -=  1
        # assert n == 0
        x = _c(s, y0, -x, _0_0)
        c = _c(c, z0, z1, _1_0)

        # Gauss-Schreiber to Gauss-Krueger TM
        # C{cmath.polar} handles INF, NAN, etc.
        k, g = polar(c)
        g = degrees(g)
    else:  # E.isSpherical
        g, k = _0_0, _1_0

    return x.real, x.imag, g, k


def _sinhcosh2(x):
    '''(INTERNAL) Like C{sincos2}.
    '''
    return sinh(x), cosh(x)


def _Xs(_Coeffs, m, E, RA=False):  # in .rhumbx
    '''(INTERNAL) Compute the C{A}, C{B} or C{RA} terms of order
       B{C{m}} for I{Krüger} series and I{rhumbx._sincosSeries},
       return a tuple with C{B{m} + 1} terms C{X}, C{X[0]==0}.
    '''
    Cs = _Coeffs[m]
    assert len(Cs) == (((m + 1) * (m + 4)) if RA else
                       ((m + 3) *  m)) // 2
    n = n_ = E.n
    if n:  # isEllipsoidal
        X = [0]  # X[0] never used, it's just an integration
        # constant, it cancels when evaluating a definite
        # integral.  Don't bother computing it, it is unused
        # in C{_Cyxgk4} above and C{rhumbx._sincosSeries}.
        _X, _p = X.append, _polynomial
        i = (m + 2) if RA else 0
        for r in _reverange(m):  # [m-1 ... 0]
            j = i + r + 1
            _X(_p(n, Cs, i, j) * n_ / Cs[j])
            i = j + 1
            n_ *= n
        X =  tuple(X)
    else:  # isSpherical
        X = _0_0s(m + 1)
    return X


# _Alp- and _BetCoeffs in .rhumbx, .rhumbBase
_AlpCoeffs = {  # Generated by Maxima on 2015-05-14 22:55:13-04:00
 4: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 4
     164, 225, -480, 360, 720,  # Alp[1]/n^1, polynomial(n), order 3
     557, -864, 390, 1440,  # Alp[2]/n^2, polynomial(n), order 2
    -1236, 427, 1680,  # PYCHOK Alp[3]/n^3, polynomial(n), order 1
     49561, 161280),  # Alp[4]/n^4, polynomial(n), order 0, count = 14
 5: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 5
    -635, 328, 450, -960, 720, 1440,  # Alp[1]/n^1, polynomial(n), order 4
     4496, 3899, -6048, 2730, 10080,  # PYCHOK Alp[2]/n^2, polynomial(n), order 3
     15061, -19776, 6832, 26880,  # PYCHOK Alp[3]/n^3, polynomial(n), order 2
    -171840, 49561, 161280,  # Alp[4]/n^4, polynomial(n), order 1
     34729, 80640),  # PYCHOK Alp[5]/n^5, polynomial(n), order 0, count = 20
 6: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 6
     31564, -66675, 34440, 47250, -100800, 75600, 151200,  # Alp[1]/n^1, polynomial(n), order 5
    -1983433, 863232, 748608, -1161216, 524160, 1935360,  # PYCHOK Alp[2]/n^2, polynomial(n), order 4
     670412, 406647, -533952, 184464, 725760,  # Alp[3]/n^3, polynomial(n), order 3
     6601661, -7732800, 2230245, 7257600,  # Alp[4]/n^4, polynomial(n), order 2
    -13675556, 3438171, 7983360,  # PYCHOK Alp[5]/n^5, polynomial(n), order 1
     212378941, 319334400),  # Alp[6]/n^6, polynomial(n), order 0, count = 27
 7: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 7
     1804025, 2020096, -4267200, 2204160, 3024000, -6451200, 4838400, 9676800,  # Alp[1]/n^1, polynomial(n), order 6
     4626384, -9917165, 4316160, 3743040, -5806080, 2620800, 9676800,  # Alp[2]/n^2, polynomial(n), order 5
    -67102379, 26816480, 16265880, -21358080, 7378560, 29030400,  # PYCHOK Alp[3]/n^3, polynomial(n), order 4
     155912000, 72618271, -85060800, 24532695, 79833600,  # Alp[4]/n^4, polynomial(n), order 3
     102508609, -109404448, 27505368, 63866880,  # Alp[5]/n^5, polynomial(n), order 2
    -12282192400, 2760926233, 4151347200,  # PYCHOK Alp[6]/n^6, polynomial(n), order 1
     1522256789, 1383782400),  # Alp[7]/n^7, polynomial(n), order 0, count = 35
 8: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 8
    -75900428, 37884525, 42422016, -89611200, 46287360, 63504000, -135475200, 101606400, 203212800,  # Alp[1]/n^1, polynomial(n), order 7
     148003883, 83274912, -178508970, 77690880, 67374720, -104509440, 47174400, 174182400,  # PYCHOK Alp[2]/n^2, polynomial(n), order 6
     318729724, -738126169, 294981280, 178924680, -234938880, 81164160, 319334400,  # PYCHOK Alp[3]/n^3, polynomial(n), order 5
    -40176129013, 14967552000, 6971354016, -8165836800, 2355138720, 7664025600,  # Alp[4]/n^4, polynomial(n), order 4
     10421654396, 3997835751, -4266773472, 1072709352, 2490808320,  # PYCHOK Alp[5]/n^5, polynomial(n), order 3
     175214326799, -171950693600, 38652967262, 58118860800,  # PYCHOK Alp[6]/n^6, polynomial(n), order 2
    -67039739596, 13700311101, 12454041600,  # PYCHOK Alp[7]/n^7, polynomial(n), order 1
     1424729850961, 743921418240)  # PYCHOK Alp[8]/n^8, polynomial(n), order 0, count = 44
}
_B1Coeffs = {  # Generated by Maxima on 2015-05-14 22:55:13-04:00
 2: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER/2 == 2
     1, 16, 64, 64),  # b1 * (n + 1), polynomial(n2), order 2
 3: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER/2 == 3
     1, 4, 64, 256, 256),  # b1 * (n + 1), polynomial(n2), order 3
 4: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER/2 == 4
     25, 64, 256, 4096, 16384, 16384)  # PYCHOK b1 * (n + 1), polynomial(n2), order 4
}
_BetCoeffs = {  # Generated by Maxima on 2015-05-14 22:55:13-04:00
 4: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 4
    -4, 555, -960, 720, 1440,  # Bet[1]/n^1, polynomial(n), order 3
    -437, 96, 30, 1440,  # Bet[2]/n^2, polynomial(n), order 2
    -148, 119, 3360,  # Bet[3]/n^3, polynomial(n), order 1
     4397, 161280),  # PYCHOK Bet[4]/n^4, polynomial(n), order 0, count = 14
 5: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 5
    -3645, -64, 8880, -15360, 11520, 23040,  # Bet[1]/n^1, polynomial(n), order 4
     4416, -3059, 672, 210, 10080,  # PYCHOK Bet[2]/n^2, polynomial(n), order 3
    -627, -592, 476, 13440,  # Bet[3]/n^3, polynomial(n), order 2
    -3520, 4397, 161280,  # Bet[4]/n^4, polynomial(n), order 1
     4583, 161280),  # PYCHOK Bet[5]/n^5, polynomial(n), order 0, count = 20
 6: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 6
     384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,  # Bet[1]/n^1, polynomial(n), order 5
    -1118711, 1695744, -1174656, 258048, 80640, 3870720,  # PYCHOK Bet[2]/n^2, polynomial(n), order 4
     22276, -16929, -15984, 12852, 362880,  # Bet[3]/n^3, polynomial(n), order 3
    -830251, -158400, 197865, 7257600,  # PYCHOK Bet[4]/n^4, polynomial(n), order 2
    -435388, 453717, 15966720,  # PYCHOK Bet[5]/n^5, polynomial(n), order 1
     20648693, 638668800),  # Bet[6]/n^6, polynomial(n), order 0, count = 27
 7: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 7
    -5406467, 6156736, -6123600, -107520, 14918400, -25804800, 19353600, 38707200,  # Bet[1]/n^1, polynomial(n), order 6
     829456, -5593555, 8478720, -5873280, 1290240, 403200, 19353600,  # PYCHOK Bet[2]/n^2, polynomial(n), order 5
     9261899, 3564160, -2708640, -2557440, 2056320, 58060800,  # PYCHOK Bet[3]/n^3, polynomial(n), order 4
     14928352, -9132761, -1742400, 2176515, 79833600,  # PYCHOK Bet[4]/n^4, polynomial(n), order 3
    -8005831, -1741552, 1814868, 63866880,  # Bet[5]/n^5, polynomial(n), order 2
    -261810608, 268433009, 8302694400,  # Bet[6]/n^6, polynomial(n), order 1
     219941297, 5535129600),  # PYCHOK Bet[7]/n^7, polynomial(n), order 0, count = 35
 8: (  # GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 8
     31777436, -37845269, 43097152, -42865200, -752640, 104428800, -180633600, 135475200, 270950400,  # Bet[1]/n^1, polynomial(n), order 7
     24749483, 14930208, -100683990, 152616960, -105719040, 23224320, 7257600, 348364800,  # Bet[2]/n^2, polynomial(n), order 6
    -232468668, 101880889, 39205760, -29795040, -28131840, 22619520, 638668800,  # PYCHOK Bet[3]/n^3, polynomial(n), order 5
     324154477, 1433121792, -876745056, -167270400, 208945440, 7664025600,  # Bet[4]/n^4, polynomial(n), order 4
     457888660, -312227409, -67920528, 70779852, 2490808320,    # Bet[5]/n^5, polynomial(n), order 3
    -19841813847, -3665348512, 3758062126, 116237721600,  # PYCHOK Bet[6]/n^6, polynomial(n), order 2
    -1989295244, 1979471673, 49816166400,  # PYCHOK Bet[7]/n^7, polynomial(n), order 1
     191773887257, 3719607091200)  # Bet[8]/n^8, polynomial(n), order 0, count = 44
}

assert set(_AlpCoeffs.keys()) == set(_BetCoeffs.keys())

if __name__ == '__main__':

    from pygeodesy.interns import _usage
    from sys import argv, exit as _exit

    _exit(_usage(*argv).replace('.ktm', '.etm -series'))

# **) MIT License
#
# Copyright (C) 2022-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
