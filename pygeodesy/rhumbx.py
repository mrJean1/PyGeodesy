
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ classes U{Rhumb
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>} from
I{GeographicLib version 2.0}.

Class L{RhumbLine} has been enhanced with methods C{intersection2} and C{nearestOn4} to iteratively
find the intersection of two rhumb lines, respectively the nearest point on a rumb line along a
geodesic or perpendicular rhumb line.

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

from pygeodesy.basics import copysign0, neg, unsigned0, _zip
from pygeodesy.constants import NAN, PI_2, _0_0s, _0_0, _0_5, \
                               _1_0, _2_0, _4_0, _720_0, _over
# from pygeodesy.ellipsoids import _EWGS84  # from .karney
from pygeodesy.errors import itemsorted, RhumbError, _Xorder
from pygeodesy.fmath import hypot, hypot1
# from pygeodesy.fsums import fsum1f_  # _MODS
from pygeodesy.interns import NN, _COMMASPACE_
from pygeodesy.karney import _atan2d, Caps, _diff182, GDict, _GTuple, \
                             _norm180,  _EWGS84
from pygeodesy.ktm import KTransverseMercator, _Xs, \
                         _AlpCoeffs, _BetCoeffs  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.props import deprecated_method, Property, Property_RO, property_RO
from pygeodesy.rhumbBase import RhumbBase, RhumbLineBase,  Int, pairs, \
                                sincos2_, _update_all_rls
# from pygeodesy.streprs import pairs  # from .rhumbBase
# from pygeodesy.units import Int  # from .rhumbBase
# from pygeodesy.utily import sincos2_  # from .rhumbBase

from math import asinh, atan, cos, cosh, fabs, radians, sin, sinh, sqrt, tan

__all__ = _ALL_LAZY.rhumbx
__version__ = '23.09.15'


class Rhumb(RhumbBase):
    '''Class to solve the I{direct} and I{inverse rhumb} problems, based on
       I{elliptic functions} or I{Krüger} series expansion.

       @see: The U{Detailed Description<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Rhumb.html>} of I{Karney}'s C++ C{Rhumb Class}.
    '''
    _mRA = 6  # see .RAorder

    def __init__(self, a_earth=_EWGS84, f=None, exact=True, name=NN, **RA_TMorder):
        '''New C{rhumbx.Rhumb}.

           @kwarg a_earth: This rhumb's earth model (L{Ellipsoid}, L{Ellipsoid2},
                           L{a_f2Tuple}, L{Datum}, 2-tuple C{(a, f)}) or the
                           (equatorial) radius (C{scalar}).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     a C{scalar}, ignored otherwise.
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

    def Direct(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Solve the I{direct rhumb} problem, optionally with the area.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees}).
           @arg s12: Distance along the rhumb line from the given to
                     the destination point (C{meter}), can be negative.

           @return: L{GDict} with 2 up to 8 items C{lat2, lon2, a12, S12,
                    lat1, lon1, azi12, s12} with the destination point's
                    latitude C{lat2} and longitude C{lon2} in C{degrees},
                    the rhumb angle C{a12} in C{degrees} and area C{S12}
                    under the rhumb line in C{meter} I{squared}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses
                  a pole, the longitude of the second point is indeterminate
                  and C{NAN} is returned for C{lon2} and area C{S12}.

           @note: If the given point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        rl = RhumbLine(self, lat1, lon1, azi12, caps=Caps.LINE_OFF,
                                                name=self.name)
        return rl.Position(s12, outmask | self._debug)  # lat2, lon2, S12

    @deprecated_method
    def Direct7(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):
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
        return self._E._elliptic_e12  # _MODS.elliptic.Elliptic(-self._E._e12)

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH_DISTANCE):
        '''Solve the I{inverse rhumb} problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).

           @return: L{GDict} with 5 to 8 items C{azi12, s12, a12, S12,
                    lat1, lon1, lat2, lon2}, the rhumb line's azimuth C{azi12}
                    in compass C{degrees} between C{-180} and C{+180}, the
                    distance C{s12} and rhumb angle C{a12} between both points
                    in C{meter} respectively C{degrees} and the area C{S12}
                    under the rhumb line in C{meter} I{squared}.

           @note: The shortest rhumb line is found.  If the end points are
                  on opposite meridians, there are two shortest rhumb lines
                  and the East-going one is chosen.

           @note: If either point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        r, Cs = GDict(name=self.name), Caps
        if (outmask & Cs.AZIMUTH_DISTANCE_AREA):
            r.set_(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2)
            E = self.ellipsoid
            psi1  = E.auxIsometric(lat1)
            psi2  = E.auxIsometric(lat2)
            psi12 = psi2 - psi1
            lon12, _ = _diff182(lon1, lon2)
            if (outmask & Cs.AZIMUTH):
                r.set_(azi12=_atan2d(lon12, psi12))
            if (outmask & Cs.DISTANCE):
                a12 = hypot(lon12, psi12) * self._DIsometric2Rectifyingd(psi2, psi1)
                s12 = a12 * E._L_90
                r.set_(s12=s12, a12=copysign0(a12, s12))
            if (outmask & Cs.AREA):
                r.set_(S12=self._S12d(lon12, psi2, psi1))
            if ((outmask | self._debug) & Cs._DEBUG_INVERSE):  # PYCHOK no cover
                r.set_(a=E.a, f=E.f, f1=E.f1, L=E.L,
                       b=E.b, e=E.e, e2=E.e2, k2=self._eF.k2,
                       lon12=lon12, psi1=psi1, exact=self.exact,
                       psi12=psi12, psi2=psi2)
        return r

#   def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
#       '''Return the distance in C{meter} and the forward and
#          reverse azimuths (initial and final bearing) in C{degrees}.
#
#          @return: L{Distance3Tuple}C{(distance, initial, final)}.
#       '''
#       r = self.Inverse(lat1, lon1, lat2, lon2)
#       return Distance3Tuple(r.s12, r.azi12, r.azi12)

    @deprecated_method
    def Inverse7(self, lat1, lon1, azi12, s12, outmask=Caps.AZIMUTH_DISTANCE_AREA):
        '''DEPRECATED, use method L{Rhumb.Inverse8}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self.Inverse8(lat1, lon1, azi12, s12, outmask=outmask)._to7Tuple()

    def _meanSinXi(self, x, y):  # radians
        s = _Dlog(cosh(x), cosh(y)) * _Dcosh(x, y)
        if self.f:
            s += _sincosSeries(False, _gd(x), _gd(y), *self._RA2) * _Dgd(x, y)
        return s

    @deprecated_method
    def orders(self, RAorder=None, TMorder=None):  # PYCHOK expected
        '''DEPRECATED, use properties C{RAorder} and/or C{TMorder}.

           Get and set the I{RAorder} and/or I{TMorder}.

           @kwarg RAorder: I{Rhumb Area} order (C{int}, 4, 5, 6, 7
                           or 8).
           @kwarg TMorder: I{Transverse Mercator} order (C{int}, 4,
                           5, 6, 7 or 8).

           @return: L{RhumbOrder2Tuple}C{(RAorder, TMorder)} with
                    the previous C{RAorder} and C{TMorder} setting.
        '''
        t = RhumbOrder2Tuple(self.RAorder, self.TMorder)
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
        n = _Xorder(_RACoeffs, RhumbError, RAorder=order)
        if self._mRA != n:
            _update_all_rls(self)
            self._mRA = n

    @Property_RO
    def _RhumbLine(self):
        '''(INTERNAL) Get this module's C{RhumbLine} class.
        '''
        return RhumbLine

    def _S12d(self, lon12, psi2, psi1):  # degrees
        '''(INTERNAL) Compute the area C{S12}.
        '''
        r = (self.ellipsoid.areax if self.exact else
             self.ellipsoid.area) * lon12 / _720_0
        r *= self._meanSinXi(radians(psi2), radians(psi1))
        return r

    @Property
    def TMorder(self):
        '''Get the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mTM

    @TMorder.setter  # PYCHOK setter!
    def TMorder(self, order):
        '''Set the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).

           @note: Setting C{TMorder} turns property C{exact} off.
        '''
        self.exact = self._TMorder(order)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{Rhumb} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: Tuple items (C{str}).
        '''
        d = dict(ellipsoid=self.ellipsoid, RAorder=self.RAorder,
                     exact=self.exact,     TMorder=self.TMorder)
        return sep.join(pairs(itemsorted(d, asorted=False), prec=prec))


class RhumbLine(RhumbLineBase):
    '''Compute one or several points on a single rhumb line.

       Class C{RhumbLine} facilitates the determination of points on
       a single rhumb line.  The starting point (C{lat1}, C{lon1})
       and the azimuth C{azi12} are specified once.

       Method C{RhumbLine.Position} returns the location of an other
       point at distance C{s12} along and the area C{S12} under the
       rhumb line.

       Method C{RhumbLine.intersection2} finds the intersection between
       two rhumb lines.

       Method C{RhumbLine.nearestOn4} computes the nearest point on and
       the distance to a rhumb line in different ways.
    '''
    _Rhumb = Rhumb  # rhumbx.Rhumb

    def __init__(self, rhumb, lat1=0, lon1=0, azi12=None, **caps_name):  # PYCHOK signature
        '''New C{rhumbx.RhumbLine}.

           @arg rhumb: The rhumb reference (C{rhumbx.Rhumb}).
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
    def _mu1(self):
        '''(INTERNAL) Get the I{rectifying auxiliary} latitude C{mu} (C{degrees}).
        '''
        return self.ellipsoid.auxRectifying(self.lat1)

    def Position(self, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Compute a point at a given distance on this rhumb line.

           @arg s12: The distance along this rhumb between its point and
                     the other point (C{meters}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: L{GDict} with 4 to 8 items C{azi12, a12, s12, S12, lat2,
                    lon2, lat1, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the point in C{degrees}, the rhumb angle C{a12}
                    in C{degrees} from the start point of and the area C{S12}
                    under this rhumb line in C{meter} I{squared}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses a
                  pole, the longitude of the second point is indeterminate and
                  C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        r, Cs = GDict(name=self.name), Caps
        if (outmask & Cs.LATITUDE_LONGITUDE_AREA):
            E, R = self.ellipsoid, self.rhumb
            a12  = s12 / E._L_90
            mu12 = self._calp * a12
            mu2, x90 = self._mu22(mu12, self._mu1)
            if x90:  # PYCHOK no cover
                lat2 = E.auxRectifying(mu2, inverse=True)
                lon2 = NAN
                if (outmask & Cs.AREA):
                    r.set_(S12=NAN)
            else:
                psi2 = self._psi1
                if self._calp:
                    lat2  = E.auxRectifying(mu2, inverse=True)
                    psi12 = R._DRectifying2Isometricd(mu2,
                                                self._mu1) * mu12
                    lon2  = psi12 * self._salp / self._calp
                    psi2 += psi12
                else:  # PYCHOK no cover
                    lat2 = self.lat1
                    lon2 = self._salp * s12 / self._r1rad
                if (outmask & Cs.AREA):
                    S12 = R._S12d(lon2, self._psi1, psi2)
                    r.set_(S12=unsigned0(S12))  # like .gx
                if (outmask & Cs.LONGITUDE):
                    if (outmask & Cs.LONG_UNROLL):
                        lon2 +=  self.lon1
                    else:
                        lon2  = _norm180(self._lon12 + lon2)
            r.set_(azi12=self.azi12, s12=s12, a12=a12)
            if (outmask & Cs.LATITUDE):
                r.set_(lat2=lat2, lat1=self.lat1)
            if (outmask & Cs.LONGITUDE):
                r.set_(lon2=lon2, lon1=self.lon1)
            if ((outmask | self._debug) & Cs._DEBUG_DIRECT_LINE):  # PYCHOK no cover
                r.set_(a=E.a, f=E.f, f1=E.f1, L=E.L, exact=R.exact,
                       b=E.b, e=E.e, e2=E.e2, k2=R._eF.k2,
                       calp=self._calp, mu1 =self._mu1,  mu12=mu12,
                       salp=self._salp, psi1=self._psi1, mu2=mu2)
        return r

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

    @Property_RO
    def _r1rad(self):  # PYCHOK no cover
        '''(INTERNAL) Get this rhumb line's parallel I{circle radius} (C{meter}).
        '''
        return radians(self.ellipsoid.circle4(self.lat1).radius)


class RhumbOrder2Tuple(_GTuple):
    '''2-Tuple C{(RAorder, TMorder)} with a I{Rhumb Area} and
       I{Transverse Mercator} order, both C{int}, DEPRECATED.
    '''
    _Names_ = (Rhumb.RAorder.name, Rhumb.TMorder.name)
    _Units_ = (      Int,                Int)


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
        r = (atan(d / r) if xy > -r else (atan(x) - atan(y))) / d
    else:
        r = _1_0 / r
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
    return (asinh(d / sqrt(x * y)) / d) if d else (_1_0 / x)


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

__all__ += _ALL_DOCS(Caps, Rhumb, RhumbLine)

if __name__ == '__main__':

    from pygeodesy.lazily import printf

    def _re(fmt, r3, x3):
        e3 = []
        for r, x in _zip(r3, x3):  # strict=True
            e = fabs(r - x) / fabs(x)
            e3.append('%.g' % (e,))
        printf((fmt % r3) + ' rel errors: ' + ', '.join(e3))

    # <https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve> version 2.0
    rhumb = Rhumb(exact=True)  # WGS84 default
    printf('# %r\n', rhumb)
    r = rhumb.Direct8(40.6, -73.8, 51, 5.5e6)  # from JFK about NE
    _re('# JFK NE lat2=%.8f, lon2=%.8f, S12=%.1f', (r.lat2, r.lon2, r.S12), (71.68889988, 0.25551982, 44095641862956.148438))
    r = rhumb.Inverse8(40.6, -73.8, 51.6, -0.5)  # JFK to LHR
    _re('# JFK-LHR azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (77.76838971, 5771083.383328, 37395209100030.367188))
    r = rhumb.Inverse8(40.6, -73.8, 35.8, 140.3)  # JFK to Tokyo Narita
    _re('# JFK-NRT azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (-92.388887981699639, 12782581.0676841792, -63760642939072.492))

# % python3 -m pygeodesy.rhumbx

# Rhumb(RAorder=6, TMorder=6, ellipsoid=Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e21=0.99330562, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594), exact=True)

# JFK NE lat2=71.68889988, lon2=0.25551982, S12=44095641862956.1 rel errors: 4e-11, 2e-08, 5e-16
# JFK-LHR azi12=77.76838971, s12=5771083.383 S12=37395209100030.4 rel errors: 3e-12, 5e-15, 0
# JFK-NRT azi12=-92.38888798, s12=12782581.068 S12=-63760642939072.5 rel errors: 2e-16, 3e-16, 0

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
