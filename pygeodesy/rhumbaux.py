# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ classes U{Rhumb
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>} from
I{GeographicLib version 2.2+} renamed as L{RhumbAux} respectively L{RhumbLineAux}.

Class L{RhumbLineAux} has been enhanced with methods C{intersection2} and C{nearestOn4} to iteratively
find the intersection of two rhumb lines, respectively the nearest point on a rumb line along a
geodesic or perpendicular rhumb line from an other point.

For more details, see the I{2.2} U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, especially the U{Class List<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>},
the background information on U{Rhumb lines<https://GeographicLib.SourceForge.io/C++/doc/rhumb.html>},
utility U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} and U{Online rhumb
line calculations<https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2023) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.

@note: Class L{AuxDST} requires package U{numpy<https://PyPI.org/project/numpy>} to be installed,
       version 1.16 or newer, needed only for I{exact} area calculations C{S12} in L{RhumbAux}
       and L{RhumbLineAux}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxAngle import AuxMu, AuxPhi,  hypot
from pygeodesy.auxilats.auxDLat import AuxDLat, _DClenshaw
# from pygeodesy.auxilats.auxDST import AuxDST  # _MODS
from pygeodesy.auxilats.auxily import _Dlam, _Dp0Dpsi, _Ufloats
from pygeodesy.basics import copysign0, _reverange, _zip,  _xkwds_get
from pygeodesy.constants import EPS_2, MANT_DIG, PI4, isinf, \
                               _0_0, _4_0, _720_0, _log2, _over
# from pygeodesy.ellipsoids import _EWGS84  # from .karney
# from pygeodesy.errors import itemsorted, _xkwds_get  # from .basics, ...
from pygeodesy.karney import Caps, _polynomial,  _EWGS84
# from pygeodesy.fmath import hypot  # from .auxilats.auxAngle
from pygeodesy.interns import NN, _COMMASPACE_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.props import Property_RO  # from .rhumbBase
from pygeodesy.rhumbBase import RhumbBase, RhumbLineBase, \
                                itemsorted, pairs, Property_RO
# from pygeodesy.streprs import pairs  # from .rhumbBase

from math import ceil as _ceil, fabs, radians

__all__ = _ALL_LAZY.rhumbaux
__version__ = '23.09.20'

# DIGITS = (sizeof(real) * 8) bits
#        = (ctypes.sizeof(ctypes.c_double(1.0)) * 8) bits
# For |n| <= 0.99, actual max for doubles is 2163.  This scales
# as DIGITS and for long doubles (GEOGRAPHICLIB_PRECISION = 3,
# DIGITS = 64), this becomes 2163 * 64 / 53 = 2612.  Round this
# up to 2^12 = 4096 and scale this by DIGITS//64 if DIGITS > 64.
#
# 64 = DIGITS for long double, 6 = 12 - _log2(64)
_Lbits = 1 << (int(_ceil(_log2(max(MANT_DIG, 64)))) + 6)


class RhumbAux(RhumbBase):
    '''Class to solve the I{direct} and I{inverse rhumb} problems, based
       on I{Auxiliary Latitudes} for accuracy near the poles.

       @note: Package U{numpy<https://PyPI.org/project/numpy>} must be
              installed, version 1.16 or later.
    '''

    def __init__(self, a_earth=_EWGS84, f=None, exact=True, name=NN, **TMorder):  # PYCHOK signature
        '''New C{rhumbaux.RhumbAux}.

           @kwarg a_earth: This rhumb's earth model (L{Ellipsoid}, L{Ellipsoid2},
                           L{a_f2Tuple}, L{Datum}, 2-tuple C{(a, f)}) or the
                           (equatorial) radius (C{scalar}).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     a C{scalar}, ignored otherwise.
           @kwarg exact: If C{True}, use the exact expressions for the I{Auxiliary
                         Latitudes}, otherwise use the I{Fourier} series expansion
                         (C{bool}), see also property C{exact}.
           @kwarg name: Optional name (C{str}).
           @kwarg TMorder: Optional keyword argument B{C{TMorder}}, see property
                           C{TMorder}.

           @raise ImportError: Package C{numpy} not found or not installed, only
                               required for area C{S12} when C{B{exact} is True}.

           @raise RhumbError: Invalid B{C{a_earth}}, B{C{f}} or B{C{RA_TMorder}}.
        '''
        RhumbBase.__init__(self, a_earth, f, exact, name)
        if TMorder:
            self.Tmorder = _xkwds_get(TMorder, TMorder=RhumbBase._mTM)

    def areaux(self, **exact):
        '''Get this ellipsoid's B{C{exact}} surface area (C{meter} I{squared}).

           @kwarg exact: Optional C{exact} (C{bool}), overriding this rhumb's
                         C{exact} setting, if C{True}, use the exact expression
                         for the authalic radius otherwise the I{Taylor} series.

           @return: The (signed?) surface area (C{meter} I{squared}).

           @raise AuxError: If C{B{exact}=False} and C{abs(flattening)} exceeds
                            property C{f_max}.

           @note: The area of a polygon encircling a pole can be found by adding
                  C{areaux / 2} to the sum of C{S12} for each side of the polygon.

           @see: U{The area of rhumb polygons<https://ArXiv.org/pdf/2303.03219.pdf>}
                 and method L{auxilats.AuxLat.AuthalicRadius2}.
        '''
        x = _xkwds_get(exact, exact=self.exact)
        a = (self._c2 * _720_0) if bool(x) is self.exact else \
            (self._auxD.AuthalicRadius2(exact=x, f_max=self.f_max) * PI4)
        return a

    @Property_RO
    def _auxD(self):
        return AuxDLat(self.ellipsoid)

    @Property_RO
    def _c2(self):  # radians makes _c2 a factor per degree
        return radians(self._auxD.AuthalicRadius2(exact=self.exact, f_max=self.f_max))

    def _DMu_DPsi(self, Phi1, Phi2, Chi1, Chi2):
        xD = self._auxD
        r  = xD.DRectifying(Phi1, Phi2) if self.exact else \
             xD.CRectifying(Chi1, Chi2)
        if r:
            r = _over(r, xD.DIsometric(Phi1, Phi2) if self.exact else
                        _Dlam(Chi1.tan, Chi2.tan))  # not Lambertian!
        return r

    def _Inverse4(self, lon12, r, outmask):
        '''(INTERNAL) See method C{RhumbBase.Inverse}.
        '''
        psi1, Chi1, Phi1 = self._psiChiPhi3(r.lat1)
        psi2, Chi2, Phi2 = self._psiChiPhi3(r.lat2)
        psi12 = psi2 - psi1  # radians
        lam12 = radians(lon12)
        if (outmask & Caps.DISTANCE):
            if isinf(psi1) or isinf(psi2):  # PYCHOK no cover
                d  = Phi2.toMu(self).toRadians
                d -= Phi1.toMu(self).toRadians
                s  = fabs(d)
            else:  # dmu/dpsi = dmu/dchi/dpsi/dchi
                s  = hypot(lam12, psi12)
                if s:
                    s *= self._DMu_DPsi(Phi1, Phi2, Chi1, Chi2)
            s *= self._rrm
            a = _over(s, self._mpd)
            r.set_(a12=copysign0(a, s), s12=s)
        return lam12, psi12, Chi1, Chi2

    def _latPhi2(self, mu):
        Mu  = AuxMu.fromDegrees(mu)
        Phi = Mu.toPhi(self)
        return Phi.toDegrees, Phi

    @Property_RO
    def _mpd(self):  # meter per degree
        return radians(self._rrm)  # == self.ellipsoid._Lpd

    def _psiChiPhi3(self, lat):
        Phi = AuxPhi.fromDegrees(lat)
        Chi = Phi.toChi(self)
        psi = Chi.toLambertianRadians
        return psi, Chi, Phi

    @Property_RO
    def _RA(self):  # get the coefficients for area calculation
        return tuple(_RAintegrate(self._auxD) if self.exact else
                     _RAseries(self._auxD))

#   _RhumbLine = RhumbLineAux  # see further below

    @Property_RO
    def _rrm(self):
        return self._auxD.RectifyingRadius(exact=self.exact)

    def _S12d(self, Chix, Chiy, lon12):  # degrees
        '''(INTERNAL) Compute the area C{S12} from C{._meanSinXi(Chix, Chiy) * .c2 * lon12}.
        '''
        pP, xD = self._RA, self._auxD

        tx, Phix = Chix.tan, Chix.toPhi(self)
        ty, Phiy = Chiy.tan, Chiy.toPhi(self)

        dD = xD.DParametric(Phix, Phiy) if self.exact else \
             xD.CParametric(Chix, Chiy)
        if dD:
            dD  = _over(dD, xD.DIsometric(Phix, Phiy) if self.exact else
                           _Dlam(tx, ty))  # not Lambertian!
            dD *= _DClenshaw(False, Phix.toBeta(self).normalized,
                                    Phiy.toBeta(self).normalized,
                                    pP, min(len(pP), xD.ALorder))  # Fsum
        dD += _Dp0Dpsi(tx, ty)
        dD *=  self._c2 * lon12
        return float(dD)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{Rhumb} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: Tuple items (C{str}).
        '''
        d = dict(ellipsoid=self.ellipsoid, exact=self.exact,
                   TMorder=self.TMorder)
        return sep.join(pairs(itemsorted(d, asorted=False), prec=prec))


class RhumbLineAux(RhumbLineBase):
    '''Compute one or several points on a single rhumb line.

       Class C{RhumbLineAux} facilitates the determination of points
       on a single rhumb line.  The starting point (C{lat1}, C{lon1})
       and the azimuth C{azi12} are specified once.

       Method C{RhumbLineAux.Position} returns the location of an
       other point and optionally the distance C{s12} along and the
       area C{S12} under the rhumb line.

       Method C{RhumbLineAux.intersection2} finds the intersection
       between two rhumb lines.

       Method C{RhumbLineAux.nearestOn4} computes the nearest point
       on and the distance to a rhumb line in different ways.
    '''
    _Rhumb = RhumbAux  # rhumbaux.RhumbAux

    def __init__(self, rhumb, lat1=0, lon1=0, azi12=None, **caps_name):  # PYCHOK signature
        '''New C{rhumbaux.RhumbLineAux}.

           @arg rhumb: The rhumb reference (C{rhumbaux.RhumbAux}).
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
    def _Chi1(self):
        return self._Phi1.toChi(self.rhumb)

    @Property_RO
    def _mu1(self):
        '''(INTERNAL) Get the I{rectifying auxiliary} latitude (C{degrees}).
        '''
        return self._Phi1.toMu(self.rhumb).toDegrees

    def _mu2lat(self, mu):
        '''(INTERNAL) Get the inverse I{rectifying auxiliary} latitude (C{degrees}).
        '''
        lat, _ = self.rhumb._latPhi2(mu)
        return lat

    @Property_RO
    def _Phi1(self):
        return AuxPhi.fromDegrees(self.lat1)

    def _Position4(self, a12, mu2, *unused):  # PYCHOK s12, mu2
        '''(INTERNAL) See method C{RhumbLineBase._Position}.
        '''
        R = self.rhumb
        lat2,  Phi2 = R._latPhi2(mu2)
        Chi2 = Phi2.toChi(R)
        Chi1 = self._Chi1
        lon2 = self._salp * a12
        if lon2:
            m = R._DMu_DPsi(self._Phi1, Phi2, Chi1, Chi2)
            lon2 = _over(lon2, m)
        return lat2, lon2, Chi1, Chi2

#   @Property_RO
#   def _psi1(self):
#       return self._Chi1.toLambertianRadians

RhumbAux._RhumbLine = RhumbLineAux  # PYCHOK see RhumbBase._RhumbLine


def _RAintegrate(auxD):
    # Compute coefficients by Fourier transform of integrand
    L   =  2
    fft = _MODS.auxilats.auxDST.AuxDST(L)
    f   =  auxD._qIntegrand
    c_  =  fft.transform(f)
    pP  =  []
    _P  =  pP.append
    # assert L < _Lbits
    while L < _Lbits:
        L  = fft.reset(L) * 2
        c_ = fft.refine(f, c_, _0_0)  # sentine[L]
        # assert len(c_) == L + 1
        pP[:], k = [], -1
        for j in range(1, L + 1):
            # Compute Fourier coefficients of integral
            p = (c_[j - 1] + c_[j]) / (_4_0 * j)
            if fabs(p) > EPS_2:
                k = -1  # run interrupted
            else:
                if k < 0:
                    k = 1  # mark as first small value
                if (j - k) >= ((j + 7) // 8):
                    # run of at least (j - 1) // 8 small values
                    return pP[:j]  # break while L loop
            _P(-p)
    return pP  # no convergence, use pP as-is


def _RAseries(auxD):
    # Series expansions in n for Fourier coeffients of the integral
    # @see: U{"Series expansions for computing rhumb areas"
    #       <https:#DOI.org/10.5281/zenodo.7685484>}.
    d  =  n = auxD._n
    i  =  0
    aL =  auxD.ALorder
    Cs = _RACoeffs[aL]
    # assert len(Cs) == (aL * (aL + 1)) // 2
    pP = []
    _P =  pP.append
    _p = _polynomial
    for m in _reverange(aL):  # order
        j  = i + m + 1
        _P(_p(n, Cs, i, j) * d)
        d *= n
        i  = j
    # assert i == len(pP)
    return pP


_f, _u = float, _Ufloats()
_RACoeffs = {  # Rhumb Area Coefficients in matrix Q
 4: _u(  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 4
    596 / _f(2025), -398 / _f(945), 22 / _f(45), -1 / _f(3),
    1543 / _f(4725), -118 / _f(315), 1 / _f(5),
    152 / _f(945), -17 / _f(315),
    5 / _f(252)),
 5: _u(  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 5
   -102614 / _f(467775), 596 / _f(2025), -398 / _f(945), 22 / _f(45),
   -1 / _f(3),
   -24562 / _f(155925), 1543 / _f(4725), -118 / _f(315), 1 / _f(5),
   -38068 / _f(155925), 152 / _f(945), -17 / _f(315),
   -752 / _f(10395), 5 / _f(252),
   -101 / _f(17325)),
 6: _u(  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 6
    138734126 / _f(638512875), -102614 / _f(467775), 596 / _f(2025),
    -398 / _f(945), 22 / _f(45), -1 / _f(3),
    17749373 / _f(425675250), -24562 / _f(155925), 1543 / _f(4725),
    -118 / _f(315), 1 / _f(5),
    1882432 / _f(8513505), -38068 / _f(155925), 152 / _f(945),
    -17 / _f(315),
    268864 / _f(2027025), -752 / _f(10395), 5 / _f(252),
    62464 / _f(2027025), -101 / _f(17325),
    11537 / _f(4054050)),
 7: _u(  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 7
    -565017322 / _f(1915538625), 138734126 / _f(638512875),
    -102614 / _f(467775), 596 / _f(2025), -398 / _f(945), 22 / _f(45),
    -1 / _f(3),
    -1969276 / _f(58046625), 17749373 / _f(425675250), -24562 / _f(155925),
    1543 / _f(4725), -118 / _f(315), 1 / _f(5),
    -58573784 / _f(638512875), 1882432 / _f(8513505), -38068 / _f(155925),
    152 / _f(945), -17 / _f(315),
    -6975184 / _f(42567525), 268864 / _f(2027025), -752 / _f(10395),
    5 / _f(252),
    -112832 / _f(1447875), 62464 / _f(2027025), -101 / _f(17325),
    -4096 / _f(289575), 11537 / _f(4054050),
    -311 / _f(525525)),
 8: _u(  # GEOGRAPHICLIB_RHUMBAREA_ORDER == 8
    188270561816 / _f(488462349375), -565017322 / _f(1915538625),
    138734126 / _f(638512875), -102614 / _f(467775), 596 / _f(2025),
    -398 / _f(945), 22 / _f(45), -1 / _f(3),
    2332829602 / _f(23260111875), -1969276 / _f(58046625),
    17749373 / _f(425675250), -24562 / _f(155925), 1543 / _f(4725),
    -118 / _f(315), 1 / _f(5),
    -41570288 / _f(930404475), -58573784 / _f(638512875),
    1882432 / _f(8513505), -38068 / _f(155925), 152 / _f(945),
    -17 / _f(315),
    1538774036 / _f(10854718875), -6975184 / _f(42567525),
    268864 / _f(2027025), -752 / _f(10395), 5 / _f(252),
    436821248 / _f(3618239625), -112832 / _f(1447875),
    62464 / _f(2027025), -101 / _f(17325),
    3059776 / _f(80405325), -4096 / _f(289575), 11537 / _f(4054050),
    4193792 / _f(723647925), -311 / _f(525525),
    1097653 / _f(1929727800))
}
del _f, _u, _Ufloats

__all__ += _ALL_DOCS(Caps, RhumbAux, RhumbLineAux)

if __name__ == '__main__':

    from pygeodesy import printf, RhumbAux  # PYCHOK RhumbAux is __main__.RhumbAux

    def _re(fmt, r3, x3):
        e3 = []
        for r, x in _zip(r3, x3):  # strict=True
            e = fabs(r - x) / fabs(x)
            e3.append('%.g' % (e,))
        printf((fmt % r3) + ' rel errors: ' + ', '.join(e3))

    # <https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve -p 9> version 2.2
    rhumb = RhumbAux(exact=True)  # WGS84 default
    printf('# %r\n', rhumb)
    r = rhumb.Direct8(40.6, -73.8, 51, 5.5e6)  # from JFK about NE
    _re('# JFK NE lat2=%.8f, lon2=%.8f, S12=%.1f', (r.lat2, r.lon2, r.S12), (71.688899882813, 0.2555198244234, 44095641862956.11))
    r = rhumb.Inverse8(40.6, -73.8, 51.6, -0.5)  # JFK to LHR
    _re('# JFK-LHR azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (77.7683897102557, 5771083.38332803, 37395209100030.39))
    r = rhumb.Inverse8(40.6, -73.8, 35.8, 140.3)  # JFK to Tokyo Narita
    _re('# JFK-NRT azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (-92.38888798169965, 12782581.067684170, -63760642939072.50))

# % python3 -m pygeodesy.rhumbaux

# RhumbAux(TMorder=6, ellipsoid=Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e21=0.99330562, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594), exact=True)

# JFK NE lat2=71.68889988, lon2=0.25551982, S12=44095641862956.1 rel errors: 4e-11, 2e-08, 5e-16
# JFK-LHR azi12=77.76838971, s12=5771083.383 S12=37395209100030.3 rel errors: 3e-12, 5e-15, 6e-16
# JFK-NRT azi12=-92.38888798, s12=12782581.068 S12=-63760642939072.5 rel errors: 2e-16, 3e-16, 6e-16
