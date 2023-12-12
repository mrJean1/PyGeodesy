
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s I{elliptic functions}, I{Krüger series expansion}, C++
classes U{Rhumb<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and
and U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>}
from I{GeographicLib version 2.0}, kept for backward compatibility.

Class L{RhumbLine} has been enhanced with methods C{Intersecant2}, C{Intersection} and C{PlumbTo} to
iteratively find the intersection of a rhumb line and a circle or an other rhumb line, respectively
a perpendicular geodesic or other rhumb line.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, especially the U{Class List<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>},
the background information on U{Rhumb lines<https://GeographicLib.SourceForge.io/C++/doc/rhumb.html>},
the utily U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} and U{Online
rhumb line calculations<https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2014-2022) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, neg
from pygeodesy.constants import PI_2, _0_0s, _0_0, _0_5, _1_0, \
                               _2_0, _4_0, _720_0, _over, _1_over
from pygeodesy.datums import _WGS84
# from pygeodesy.deprecated import RhumbOrder2Tuple  # _MODS
from pygeodesy.errors import RhumbError, _Xorder
from pygeodesy.fmath import hypot, hypot1
# from pygeodesy.fsums import fsum1f_  # _MODS
# from pygeodesy.interns import NN  # from .farney
from pygeodesy.karney import Caps,  NN
from pygeodesy.ktm import KTransverseMercator, _Xs, \
                         _AlpCoeffs, _BetCoeffs  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.props import deprecated_method, Property, Property_RO, property_RO
from pygeodesy.rhumb.bases import RhumbBase, RhumbLineBase,  _update_all_rls
from pygeodesy.utily import atan1, sincos2_

from math import asinh, atan, cos, cosh, radians, sin, sinh, sqrt, tan

__all__ = _ALL_LAZY.rhumb_ekx
__version__ = '23.12.03'


class Rhumb(RhumbBase):
    '''Class to solve the I{direct} and I{inverse rhumb} problems, based on
       I{elliptic functions} or I{Krüger series expansion}

       @see: The U{Detailed Description<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Rhumb.html>} of I{Karney}'s C++ C{Rhumb Class}.
    '''
    _mRA = 6  # see .RAorder

    def __init__(self, a_earth=_WGS84, f=None, exact=True, name=NN, **RA_TMorder):
        '''New C{Rhumb}.

           @kwarg a_earth: This rhumb's earth model (L{Datum}, L{Ellipsoid},
                           L{Ellipsoid2}, L{a_f2Tuple}, 2-tuple C{(a, f)}) or
                           the (equatorial) radius (C{meter}, conventionally).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     C{scalar}, ignored otherwise.
           @kwarg exact: If C{True}, use an addition theorem for elliptic integrals
                         to compute I{Divided differences}, otherwise use the I{Krüger}
                         series expansion (C{bool} or C{None}), see also properties
                         C{exact} and C{TMorder}.
           @kwarg name: Optional name (C{str}).
           @kwarg RA_TMorder: Optional keyword arguments B{C{RAorder}} and B{C{TMorder}}
                              to set the respective C{order}, see properties C{RAorder}
                              and C{TMorder} and method C{orders}.

           @raise RhumbError: Invalid B{C{a_earth}}, B{C{f}} or B{C{RA_TMorder}}.
        '''
        RhumbBase.__init__(self, a_earth, f, exact, name)
        if RA_TMorder:
            self.orders(**RA_TMorder)

    @Property_RO
    def _A2(self):  # Conformal2RectifyingCoeffs
        m = self.TMorder
        return _Xs(_AlpCoeffs, m, self.ellipsoid), m

    @Property_RO
    def _B2(self):  # Rectifying2ConformalCoeffs
        m = self.TMorder
        return _Xs(_BetCoeffs, m, self.ellipsoid), m

    def _DConformal2Rectifying(self, x, y):  # radians
        return _1_0 + (_sincosSeries(True, x, y, *self._A2) if self.f else _0_0)

    @deprecated_method
    def Direct7(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):  # PYCHOK no cover
        '''DEPRECATED, use method L{Rhumb.Direct8}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self.Direct8(lat1, lon1, azi12, s12, outmask=outmask)._to7Tuple()

    def _DIsometrict(self, phix, phiy, tphix, tphiy, _Dtan_phix_phiy):
        E = self.ellipsoid
        return _Dtan_phix_phiy   * _Dasinh(tphix, tphiy) - \
               _Dsin(phix, phiy) * _DeatanhE(sin(phix), sin(phiy), E)

    def _DIsometric2Rectifyingd(self, psix, psiy):  # degrees
        if self.exact:
            E = self.ellipsoid
            phix, phiy, tphix, tphiy = _Eaux4(E.auxIsometric, psix, psiy)
            t = _Dtant(phix - phiy, tphix, tphiy)
            r = _over(self._DRectifyingt(           tphix, tphiy, t),
                      self._DIsometrict(phix, phiy, tphix, tphiy, t))
        else:
            x, y = radians(psix), radians(psiy)
            r = self._DConformal2Rectifying(_gd(x), _gd(y)) * _Dgd(x, y)
        return r

    def _DRectifyingt(self, tphix, tphiy, _Dtan_phix_phiy):
        E = self.ellipsoid
        tbetx = E.f1 *  tphix
        tbety = E.f1 *  tphiy
        return (E.f1 * _Dtan_phix_phiy * E.b * PI_2
                     * _DfEt( tbetx, tbety, self._eF)
                     * _Datan(tbetx, tbety)) / E.L

    def _DRectifying2Conformal(self, x, y):  # radians
        return _1_0 - (_sincosSeries(True, x, y, *self._B2) if self.f else _0_0)

    def _DRectifying2Isometricd(self, mux, muy):  # degrees
        E = self.ellipsoid
        phix, phiy, tphix, tphiy = _Eaux4(E.auxRectifying, mux, muy)
        if self.exact:
            t = _Dtant(phix - phiy, tphix, tphiy)
            r = _over(self._DIsometrict(phix, phiy, tphix, tphiy, t),
                      self._DRectifyingt(           tphix, tphiy, t))
        else:
            r = self._DRectifying2Conformal(radians(mux), radians(muy)) * \
                     _Dgdinv(E.es_taupf(tphix), E.es_taupf(tphiy))
        return r

    @Property_RO
    def _eF(self):
        '''(INTERNAL) Get the ellipsoid's elliptic function.
        '''
        # .k2 = 0.006739496742276434
        return self.ellipsoid._elliptic_e12  # _MODS.elliptic.Elliptic(-self.ellipsoid._e12)

    def _Inverse4(self, lon12, r, outmask):
        '''(INTERNAL) See method C{RhumbBase.Inverse}.
        '''
        E, Cs = self.ellipsoid, Caps
        psi1  = E.auxIsometric(r.lat1)
        psi2  = E.auxIsometric(r.lat2)
        psi12 = psi2 - psi1  # degrees
        if (outmask & Cs.DISTANCE):
            a = s = hypot(lon12, psi12)
            if a:
                a *= self._DIsometric2Rectifyingd(psi2, psi1)
                s  = self._mpd * a  # == E._Lpd
                a  = copysign0(a, s)
            r.set_(a12=a, s12=s)

        if ((outmask | self._debug) & Cs._DEBUG_INVERSE):  # PYCHOK no cover
            r.set_(a=E.a, f=E.f, f1=E.f1, L=E.L,
                   b=E.b, e=E.e, e2=E.e2, k2=self._eF.k2,
                   lon12=lon12, psi1=psi1, exact=self.exact,
                   psi12=psi12, psi2=psi2)
        return lon12, psi12, psi1, psi2

    @deprecated_method
    def Inverse7(self, lat1, lon1, azi12, s12, outmask=Caps.AZIMUTH_DISTANCE_AREA):  # PYCHOK no cover
        '''DEPRECATED, use method L{Rhumb.Inverse8}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self.Inverse8(lat1, lon1, azi12, s12, outmask=outmask)._to7Tuple()

    @Property_RO
    def _mpd(self):  # meter per degree
        return self.ellipsoid._Lpd

    @Property_RO
    def _mpr(self):  # meter per radian
        return self.ellipsoid._Lpr  # degrees(._Lpd)

    @deprecated_method
    def orders(self, RAorder=None, TMorder=None):  # PYCHOK no cover
        '''DEPRECATED, use properties C{RAorder} and/or C{TMorder}.

           Get and set the I{RAorder} and/or I{TMorder}.

           @kwarg RAorder: I{Rhumb Area} order (C{int}, 4, 5, 6, 7
                           or 8).
           @kwarg TMorder: I{Transverse Mercator} order (C{int}, 4,
                           5, 6, 7 or 8).

           @return: DEPRECATED L{RhumbOrder2Tuple}C{(RAorder, TMorder)}
                    with the previous C{RAorder} and C{TMorder} setting.
        '''
        t = _MODS.deprecated.classes.RhumbOrder2Tuple(self.RAorder, self.TMorder)
        if RAorder not in (None, t.RAorder):  # PYCHOK attr
            self.RAorder = RAorder
        if TMorder not in (None, t.TMorder):  # PYCHOK attr
            self.TMorder = TMorder
        return t

    @Property_RO
    def _RA2(self):
        # for WGS84: (0, -0.0005583633519275459, -3.743803759172812e-07,  -4.633682270824446e-10,
        # RAorder 6:     -7.709197397676237e-13, -1.5323287106694307e-15, -3.462875359099873e-18)
        m = self.RAorder
        return _Xs(_RACoeffs, m, self.ellipsoid, RA=True), m

    @Property
    def RAorder(self):
        '''Get the I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mRA

    @RAorder.setter  # PYCHOK setter!
    def RAorder(self, order):
        '''Set the I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        m = _Xorder(_RACoeffs, RhumbError, RAorder=order)
        if self._mRA != m:
            _update_all_rls(self)
            self._mRA = m

#   _RhumbLine = RhumbLine  # see further below

    def _S12d(self, psi1, psi2, lon12):  # degrees
        '''(INTERNAL) Compute the area C{S12}.
        '''
        S = (self.ellipsoid.areax if self.exact else
             self.ellipsoid.area) * lon12 / _720_0
        if S:
            x, y = radians(psi1), radians(psi2)  # _meanSinXi(x, y)
            s = _Dlog(cosh(x), cosh(y)) * _Dcosh(x, y)
            if self.f:
                s += _sincosSeries(False, _gd(x), _gd(y), *self._RA2) * _Dgd(x, y)
            S *= s
        return S


class RhumbLine(RhumbLineBase):
    '''Compute one or several points on a single rhumb line.

       Class C{RhumbLine} facilitates the determination of points on
       a single rhumb line.  The starting point (C{lat1}, C{lon1})
       and the azimuth C{azi12} are specified once.
    '''
    _Rhumb = Rhumb  # rhumb.ekx.Rhumb

    def __init__(self, rhumb, lat1=0, lon1=0, azi12=None, **caps_name):  # PYCHOK signature
        '''New C{RhumbLine}.

           @arg rhumb: The rhumb reference (L{Rhumb}).
           @kwarg lat1: Latitude of the start point (C{degrees90}).
           @kwarg lon1: Longitude of the start point (C{degrees180}).
           @kwarg azi12: Azimuth of this rhumb line (compass C{degrees}).
           @kwarg caps_name: Optional keyword arguments C{B{name}=NN} and
                       C{B{caps}=0}, a bit-or'ed combination of L{Caps}
                       values specifying the required capabilities.  Include
                       C{Caps.LINE_OFF} if updates to the B{C{rhumb}} should
                       I{not} be reflected in this rhumb line.
        '''
        RhumbLineBase.__init__(self, rhumb, lat1, lon1, azi12, **caps_name)

    @Property_RO
    def _dpm12(self):  # PYCHOK no cover
        '''(INTERNAL) Get this rhumb line's parallel I{circle radius} (C{degree per meter}).
        '''
        r = self._salp
        if r:
            r = _over(r, self.ellipsoid.circle4(self.lat1).radius)
        return r

    @Property_RO
    def _mu1(self):
        '''(INTERNAL) Get the I{rectifying auxiliary} latitude (C{degrees}).
        '''
        return self.ellipsoid.auxRectifying(self.lat1)

    def _mu2lat(self, mu):
        '''(INTERNAL) Get the inverse I{rectifying auxiliary} latitude (C{degrees}).
        '''
        return self.ellipsoid.auxRectifying(mu, inverse=True)

    def _Position4(self, unused, mu2, s12, mu12):
        '''(INTERNAL) See method C{RhumbLineBase._Position}.
        '''
        psi1 = psi2 = self._psi1
        if mu12:  # self._calp != 0
            lat2  = self._mu2lat(mu2)
            psi12 = self.rhumb._DRectifying2Isometricd(mu2, self._mu1) * mu12
            lon2  = self._talp * psi12
            psi2 += psi12
        else:  # meridional
            lat2  = self.lat1
            lon2  = self._dpm12 * s12
        return lat2, lon2, psi1, psi2

    @Property_RO
    def _psi1(self):
        '''(INTERNAL) Get the I{isometric auxiliary} latitude C{psi} (C{degrees}).
        '''
        return self.ellipsoid.auxIsometric(self.lat1)

    @property_RO
    def RAorder(self):
        '''Get this rhumb line's I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self.rhumb.RAorder

Rhumb._RhumbLine = RhumbLine  # PYCHOK see RhumbBase._RhumbLine


# Use I{Divided Differences} to determine (mu2 - mu1) / (psi2 - psi1) accurately.
# Definition: _Df(x,y,d) = (f(x) - f(y)) / (x - y), @see W. M. Kahan & R. J.
# Fateman, "Symbolic computation of Divided Differences", SIGSAM Bull. 33(3),
# 7-28 (1999). U{ACM<https://DL.ACM.org/doi/pdf/10.1145/334714.334716> and @see
# U{UCB<https://www.CS.Berkeley.edu/~fateman/papers/divdiff.pdf>}, Dec 8, 1999.

def _Dasinh(x, y):
    hx = hypot1(x)
    d  = x - y
    if d:
        hx *= y
        hy  = x * hypot1(y)
        t = (d * (x + y) / (hy + hx)) if (x * y) > 0 else (hy - hx)
        r =  asinh(t) / d
    else:
        r = _1_0 / hx
    return r


def _Datan(x, y):
    xy = x  *  y
    r  = xy + _1_0
    d  = x  -  y
    if d:  # 2 * xy > -1 == 2 * xy + 1 > 0 == xy + r > 0 == xy > -r
        r = (atan1(d, r) if xy > -r else (atan1(x) - atan1(y))) / d
    else:
        r = _1_over(r)
    return r


def _Dcosh(x, y):
    return _Dsincos(x, y, sinh, sinh)


def _DeatanhE(x, y, E):  # see .albers._Datanhee
    # Deatanhe(x, y) = eatanhe((x - y) / (1 - e^2 * x * y)) / (x - y)
    e = _1_0 - E.e2 * x * y
    if e:  # assert not isnear0(e)
        d =  x - y
        e = (E._es_atanh(d / e) / d) if d else (E.e2 / e)
    return e


def _DfEt(tx, ty, eF):  # tangents
    # eF = Elliptic(-E.e12)  # -E.e2 / (1 - E.e2)
    r,  x,  y, = _1_0, atan(tx), atan(ty)
    d = x - y
    if (x * y) > 0:
        # See U{DLMF<https://DLMF.NIST.gov/19.11>}: 19.11.2 and 19.11.4
        # letting theta -> x, phi -> -y, psi -> z
        # (E(x) - E(y)) / d = E(z)/d - k2 * sin(x) * sin(y) * sin(z)/d
        # tan(z/2) = (sin(x)*Delta(y) - sin(y)*Delta(x)) / (cos(x) + cos(y))
        #          = d * Dsin(x,y) * (sin(x) + sin(y))/(cos(x) + cos(y)) /
        #                (sin(x)*Delta(y) + sin(y)*Delta(x))
        #          = t = d * Dt
        # sin(z) = 2*t/(1+t^2); cos(z) = (1-t^2)/(1+t^2)
        # Alt (this only works for |z| <= pi/2 -- however, this conditions
        # holds if x*y > 0):
        #   sin(z) = d * Dsin(x,y) * (sin(x) + sin(y)) /
        #                (sin(x)*cos(y)*Delta(y) + sin(y)*cos(x)*Delta(x))
        #   cos(z) = sqrt((1-sin(z))*(1+sin(z)))
        sx, cx, sy, cy = sincos2_(x, y)
        D  = (cx + cy) * (eF.fDelta(sy, cy) * sx +
                          eF.fDelta(sx, cx) * sy)
        D  = (sx + sy) * _Dsin(x, y) / D
        t  =  D * d
        t2 = _1_0 + t**2
        D *= _2_0 / t2
        s  =  D * d
        if s:
            c = (t + _1_0) * (_1_0 - t) / t2
            r = eF.fE(s, c, eF.fDelta(s, c)) / s
        r = D * (r - eF.k2 * sx * sy)
    elif d:
        r = (eF.fE(x) - eF.fE(y)) / d
    return r


def _Dgd(x, y):
    return _Datan(sinh(x), sinh(y)) * _Dsinh(x, y)


def _Dgdinv(x, y):  # x, y are tangents
    return _Dasinh(x, y) / _Datan(x, y)


def _Dlog(x, y):
    d = (x - y) * _0_5
    # Changed atanh(t / (x + y)) to asinh(t / (2 * sqrt(x*y))) to
    # avoid taking atanh(1) when x is large and y is 1.  This also
    # fixes bogus results being returned for the area when an endpoint
    # is at a pole.  N.B. this routine is invoked with positive x
    # and y, so the sqrt is always taken of a positive quantity.
    return (asinh(d / sqrt(x * y)) / d) if d else _1_over(x)


def _Dsin(x, y):
    return _Dsincos(x, y, sin, cos)


def _Dsincos(x, y, sin_, cos_):
    r = cos_((x + y) * _0_5)
    d =      (x - y) * _0_5
    if d:
        r *= sin_(d) / d
    return r


def _Dsinh(x, y):
    return _Dsincos(x, y, sinh, cosh)


def _Dtan(x, y):  # PYCHOK no cover
    return _Dtant(x - y, tan(x), tan(y))


def _Dtant(dxy, tx, ty):
    txy = tx  * ty
    r   = txy + _1_0
    if dxy:  # 2 * txy > -1 == 2 * txy + 1 > 0 == txy + r > 0 == txy > -r
        r = ((tan(dxy) * r) if txy > -r else (tx - ty)) / dxy
    return r


def _Eaux4(E_aux, mu_psi_x, mu_psi_y):  # degrees
    # get inverse auxiliary lats in radians and tangents
    phix = radians(E_aux(mu_psi_x, inverse=True))
    phiy = radians(E_aux(mu_psi_y, inverse=True))
    return phix, phiy, tan(phix), tan(phiy)


def _gd(x):
    return atan(sinh(x))


def _sincosSeries(sinp, x, y, C, n):
    # N.B. C[] has n+1 elements of which
    #      C[0] is ignored and n >= 0
    # Use Clenshaw summation to evaluate
    #   m = (g(x) + g(y)) / 2        -- mean value
    #   s = (g(x) - g(y)) / (x - y)  -- average slope
    # where
    #   g(x) = sum(C[j] * SC(2 * j * x), j = 1..n)
    #   SC = sinp ? sin : cos
    #   CS = sinp ? cos : sin
    # ...
    d, _neg = (x - y), neg
    sp, cp, sd, cd = sincos2_(x + y, d)
    sd = (sd / d) if d else _1_0
    s  = _neg(sp * sd)  # negative
    # 2x2 matrices in row-major order
    a1 = s  *  d**2
    a2 = s  * _4_0
    a0 = a3 = _2_0 * cp * cd  # m
    b2 = b1 = _0_0s(4)
    if n > 0:
        b1 = C[n], _0_0, _0_0, C[n]

    _fsum = _MODS.fsums.fsum1f_
    for j in range(n - 1, 0, -1):  # C[0] unused
        b1, b2, Cj = b2, b1, C[j]
        # b1 = a * b2 - b1 + C[j] * I
        m0, m1, m2, m3 = b2
        n0, n1, n2, n3 = map(_neg, b1)
        b1 = (_fsum(a0 * m0, a1 * m2, n0, Cj),
              _fsum(a0 * m1, a1 * m3, n1),
              _fsum(a2 * m0, a3 * m2, n2),
              _fsum(a2 * m1, a3 * m3, n3, Cj))
    # Here are the full expressions for m and s
    # f01, f02, f11, f12 = (0, 0, cd * sp,  2 * sd * cp) if sinp else \
    #                      (1, 0, cd * cp, -2 * sd * sp)
    # m = -b2[1] * f02 + (C[0] - b2[0]) * f01 + b1[0] * f11 + b1[1] * f12
    # s = -b2[2] * f01 + (C[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12
    cd *=  b1[2]
    sd *=  b1[3] * _2_0
    s   = _fsum(cd * sp,      sd * cp) if sinp else \
          _fsum(cd * cp, _neg(sd * sp), _neg(b2[2]))
    return s


_RACoeffs = {  # Generated by Maxima on 2015-05-15 08:24:04-04:00
 4: (  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 4
     691, 7860, -20160, 18900, 0, 56700,  # R[0]/n^0, polynomial(n), order 4
     1772, -5340, 6930, -4725, 14175,  # R[1]/n^1, polynomial(n), order 3
    -1747, 1590, -630, 4725,  # PYCHOK R[2]/n^2, polynomial(n), order 2
     104, -31, 315,  # R[3]/n^3, polynomial(n), order 1
    -41, 420),  # PYCHOK R[4]/n^4, polynomial(n), order 0, count = 20
 5: (  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 5
    -79036, 22803, 259380, -665280, 623700, 0, 1871100,  # PYCHOK R[0]/n^0, polynomial(n), order 5
     41662, 58476, -176220, 228690, -155925, 467775,  # PYCHOK R[1]/n^1, polynomial(n), order 4
     18118, -57651, 52470, -20790, 155925,  # PYCHOK R[2]/n^2, polynomial(n), order 3
    -23011, 17160, -5115, 51975,  # PYCHOK R[3]/n^3, polynomial(n), order 2
     5480, -1353, 13860,  # PYCHOK R[4]/n^4, polynomial(n), order 1
    -668, 5775),    # PYCHOK R[5]/n^5, polynomial(n), order 0, count = 27
 6: (  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 6
     128346268, -107884140, 31126095, 354053700, -908107200, 851350500, 0, 2554051500,  # R[0]/n^0, polynomial(n), order 6
    -114456994, 56868630, 79819740, -240540300, 312161850, -212837625, 638512875,  # PYCHOK R[1]/n^1, polynomial(n), order 5
     51304574, 24731070, -78693615, 71621550, -28378350, 212837625,  # R[2]/n^2, polynomial(n), order 4
     1554472, -6282003, 4684680, -1396395, 14189175,  # R[3]/n^3, polynomial(n), order 3
    -4913956, 3205800, -791505, 8108100,  # PYCHOK R[4]/n^4, polynomial(n), order 2
     1092376, -234468, 2027025,  # R[5]/n^5, polynomial(n), order 1
    -313076, 2027025),    # PYCHOK R[6]/n^6, polynomial(n), order 0, count = 35
 7: (  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 7
    -317195588, 385038804, -323652420, 93378285, 1062161100, -2724321600, 2554051500, 0, 7662154500,  # PYCHOK R[0]/n^0, polynomial(n), order 7
     258618446, -343370982, 170605890, 239459220, -721620900, 936485550, -638512875, 1915538625,  # PYCHOK R[1]/n^1, polynomial(n), order 6
    -248174686, 153913722, 74193210, -236080845, 214864650, -85135050, 638512875,  # PYCHOK R[2]/n^2, polynomial(n), order 5
     114450437, 23317080, -94230045, 70270200, -20945925, 212837625,  # PYCHOK R[3]/n^3, polynomial(n), order 4
     15445736, -103193076, 67321800, -16621605, 170270100,  # PYCHOK R[4]/n^4, polynomial(n), order 3
    -27766753, 16385640, -3517020, 30405375,  # PYCHOK R[4]/n^4, polynomial(n), order 3
     4892722, -939228, 6081075,  # PYCHOK R[4]/n^4, polynomial(n), order 3
    -3189007, 14189175),    # PYCHOK R[7]/n^7, polynomial(n), order 0, count = 44
 8: (  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 8
     71374704821, -161769749880, 196369790040, -165062734200, 47622925350, 541702161000, -1389404016000, 1302566265000, 0, 3907698795000,  # R[0]/n^0, polynomial(n), order 8
    -13691187484, 65947703730, -87559600410, 43504501950, 61062101100, -184013329500, 238803815250, -162820783125, 488462349375,  # PYCHOK R[1]/n^1, polynomial(n), order 7
     30802104839, -63284544930, 39247999110, 18919268550, -60200615475, 54790485750, -21709437750, 162820783125,  # R[2]/n^2, polynomial(n), order 6
    -8934064508, 5836972287, 1189171080, -4805732295, 3583780200, -1068242175, 10854718875,  # PYCHOK R[3]/n^3, polynomial(n), order 5
     50072287748, 3938662680, -26314234380, 17167059000, -4238509275, 43418875500,  # R[4]/n^4, polynomial(n), order 4
     359094172, -9912730821, 5849673480, -1255576140, 10854718875,  # R[5]/n^5, polynomial(n), order 3
    -16053944387, 8733508770, -1676521980, 10854718875,  # PYCHOK R[6]/n^6, polynomial(n), order 2
     930092876, -162639357, 723647925,  # R[7]/n^7, polynomial(n), order 1
    -673429061, 1929727800)    # PYCHOK R[8]/n^8, polynomial(n), order 0, count = 54
}

__all__ += _ALL_DOCS(Caps)

# **) MIT License
#
# Copyright (C) 2022-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
