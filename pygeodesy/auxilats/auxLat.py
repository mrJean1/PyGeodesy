
# -*- coding: utf-8 -*-

u'''Class L{AuxLat} transcoded to Python from I{Karney}'s C++ class U{AuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2023) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.

@see: U{Auxiliary latitudes<https:#GeographicLib.SourceForge.io/C++/doc/auxlat.html>}
      U{On auxiliary latitudes<https:#ArXiv.org/abs/2212.05818>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxAngle import AuxAngle, AuxBeta, AuxChi, _AuxClass, \
                                        AuxMu, AuxPhi, AuxTheta, AuxXi
from pygeodesy.auxilats.auxily import Aux, _sc, _sn, _Ufloats
from pygeodesy.basics import isscalar, _reverange, _xinstanceof
from pygeodesy.constants import INF, MAX_EXP, MIN_EXP, NAN, PI_2, PI_4, _EPSqrt, \
                               _0_0, _0_0s, _0_1, _0_25, _0_5, _1_0, _2_0, _3_0, \
                               _360_0, isfinite, isinf, isnan, _log2, _over
from pygeodesy.datums import _ellipsoidal_datum, _WGS84,  Ellipsoid
# from pygeodesy.ellipsoids import Ellipsoid  # from .datums
from pygeodesy.elliptic import Elliptic as _Ef
from pygeodesy.errors import AuxError, _xkwds, _xkwds_get, _Xorder
# from pygeodesy.fmath import cbrt  # from .karney
from pygeodesy.fsums import Fsum, _sum
from pygeodesy.karney import _2cos2x,  _ALL_DOCS, cbrt, _MODS
from pygeodesy.interns import NN, _DOT_, _UNDER_  # _earth_
# from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS  # from .karney
from pygeodesy.props import Property, Property_RO, _update_all
from pygeodesy.units import Degrees, Meter
# from pygeodesy.utily import _passarg  # _MODS

from math import asinh, atan, atan2, copysign, cosh, fabs, sin, sinh, sqrt
try:
    from math import exp2 as _exp2
except ImportError:  # Python 3.11-

    def _exp2(x):
        return pow(_2_0, x)

__all__ = ()
__version__ = '23.09.14'

_TRIPS = 1024  # XXX 2 or 3?


class AuxLat(AuxAngle):
    '''Accurate conversion between I{Auxiliary} latitudes on an ellipsoid.

       Latitudes are represented by instances of L{AuxAngle} in order to
       maintain precision close to the poles, I{Authalic} latitude C{Xi},
       I{Conformal} C{Chi}, I{Geocentric} C{Theta}, I{Geographic} C{Phi},
       I{Parametric} C{Beta} and I{Rectifying} C{Mu}.

       @see: I{Karney}'s C++ class U{AuxLatitude
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}.
    '''
    _csc =  dict()  # global coeffs cache: [aL][k], upto max(k) * (4 + 6 + 8) floats
    _E   = _WGS84.ellipsoid
    _mAL =  6       # 4, 6 or 8 aka Lmax

    def __init__(self, a_earth=_WGS84, f=None, b=None, name=NN, **ALorder):
        '''New L{AuxLat} instance on an ellipsoid or datum.

           @arg a_earth: Equatorial radius, semi-axis (C{meter}) or an
                         ellipsoid or datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg f: Flattening: M{(a - b) / a} (C{float}, near zero for
                     spherical), ignored if B{C{a_earth}} is not scalar.
           @kwarg b: Optional polar radius, semi-axis (C{meter}, same
                     units as B{C{a_earth}}), ignored if B{C{a_earth}}
                     is not scalar.
           @kwarg ALorder: Optional keyword arguments B{C{ALorder}} to
                           set L{AuxLat}'s C{order}, see property
                           C{ALorder}.
           @kwarg name: Optional name (C{str}).
        '''
        if a_earth is not _WGS84:
            n = name or AuxLat.__name__
            try:
                if b is f is None:
                    E = _ellipsoidal_datum(a_earth, name=n).ellipsoid  # XXX raiser=_earth_
                elif isscalar(a_earth):
                    E =  Ellipsoid(a_earth, f=f, b=b, name=_UNDER_(n))
                else:
                    raise ValueError()
            except Exception as x:
                raise AuxError(a_earth=a_earth, f=f, b=b, cause=x)
            self._E = E

        if name:
            self.name = name
        if ALorder:
            self.ALorder = _xkwds_get(ALorder, ALorder=AuxLat._mAL)

    @Property_RO
    def a(self):
        '''Get the C{ellipsoid}'s equatorial radius (C{meter}, conventionally).
        '''
        return self.ellipsoid.a

    equatoradius = a

    @Property
    def ALorder(self):
        '''Get the I{AuxLat} order (C{int}, 4, 6 or 8).
        '''
        return self._mAL

    @ALorder.setter  # PYCHOK setter!
    def ALorder(self, order):
        '''Set the I{AuxLat} order (C{int}, 4, 6 or 8).
        '''
        m = _Xorder(_AR2Coeffs, AuxError, ALorder=order)
        if self._mAL != m:
            _update_all(self)
            self._mAL = m

    def _atanhee(self, tphi):  # see Ellipsoid._es_atanh, .albers._atanhee
        # atanh(e * sphi) = asinh(e' * sbeta)
        f =  self.f
        s = _sn(self._fm1 * tphi) if f > 0 else _sn(tphi)
        if f:  # atanh(e * sphi) = asinh(e' * sbeta)
            e = self._e
            s = (atan(e * s) if f < 0 else asinh(self._e1 * s)) / e
        return s

    def Authalic(self, Phi, **diff_name):
        '''Convert I{Geographic} to I{Aunthalic} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Parametric latitude, C{Beta} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        tphi = fabs(Phi.tan)
        if isfinite(tphi) and tphi and self.f:
            y, x = Phi._yx_normalized
            q    = self._q
            qv   = self._qf(tphi)
            Dq2  = self._Dq(tphi)
            Dq2 *= (q + qv) / (fabs(y) + _1_0)  # _Dq(-tphi)
            Xi   = AuxXi(copysign(qv, Phi.y), x * sqrt(Dq2),
                         name=_xkwds_get(diff_name, name=Phi.name))

            if _xkwds_get(diff_name, diff=False):
                if isnan(tphi):
                    d = self._e2m1_sq2
                else:
                    c  =  self.Parametric(Phi)._x_normalized
                    d  = _over(c, Xi._x_normalized)**3
                    d *= _over(c, x) * _over(_2_0, q)
                Xi._diff = d
        else:
            Xi = AuxXi(*Phi._yx)  # diff default
        # assert Xi._AUX == Aux.XI
        return Xi

    def AuthalicRadius2(self, exact=False, f_max=_0_1):
        '''Get the I{Authalic} radius I{squared}.

           @kwarg exact: If C{True}, use the exact expression, otherwise
                         the I{Taylor} series.
           @kwarg f_max: C{Flattening} not to exceed (C{float}).

           @return: Authalic radius I{squared} (C{meter} I{squared}, same
                    units as the ellipsoid axes).

           @raise AuxError: If C{B{exact}=False} and C{abs(flattening)}
                            exceeds C{f_max}.
        '''
        f = self.f
        if exact or not f:
            c2 = self.ellipsoid.b2 * self._q  # == ellipsoid.c2x * 2
        elif fabs(f) < f_max:
            # Using a * (a + b) / 2 as the multiplying factor leads to a rapidly
            # converging series in n.  Of course, using this series isn't really
            # necessary, since the exact expression is simple to evaluate.  However,
            # we do it for consistency with RectifyingRadius and, presumably, the
            # roundoff error is smaller compared to that for the exact expression.
            m   =  self.ALorder
            c2  = _MODS.karney._polynomial(self._n, _AR2Coeffs[m], 0, m)
            c2 *=  self.a * (self.a + self.b)
        else:
            raise AuxError(exact=exact, f=f, f_max=f_max)
        return c2 * _0_5

    @Property_RO
    def b(self):
        '''Get the C{ellipsoid}'s polar radius (C{meter}, conventionally).
        '''
        return self.ellipsoid.b  # a * (_1_0 - f)

    polaradius = b

    def _coeffs(self, auxout, auxin):
        # Get the polynomial coefficients as 4-, 6- or 8-tuple
        aL = self.ALorder  # aka Lmax
        k  = Aux._1d(auxout, auxin)
        try:
            cs = AuxLat._csc[aL][k]
        except KeyError:
            if auxout == auxin:
                cs = _0_0s(aL)  # not cached
            else:
                Cx = _CXcoeffs(aL)
                try:
                    Cx = Cx[auxout][auxin]
                except KeyError as x:
                    raise AuxError(auxout=auxout, auxin=auxin, cause=x)

                d = x = n = self._n
                if Aux.use_n2(auxin) and Aux.use_n2(auxout):
                    x = self._n2

                    def _m(m):
                        return m // 2
                else:
                    def _m(m):  # PYCHOK expected
                        return m

                i  = 0
                cs = []
                _p = _MODS.karney._polynomial
                for r in _reverange(aL):
                    j  = i + _m(r) + 1  # order m = j - i - 1
                    cs.append(_p(x, Cx, i, j) * d)
                    d *= n
                    i  = j
                # assert i == len(Cx) and len(cs) == aL
                AuxLat._csc.setdefault(aL, {})[k] = tuple(cs)
        return cs

    def Conformal(self, Phi, **diff_name):
        '''Convert I{Geographic} to I{Conformal} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Conformal latitude, C{Chi} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        tphi = tchi = fabs(Phi.tan)
        if isfinite(tphi) and tphi and self.f:
            sig   =  sinh(self._e2 * self._atanhee(tphi))
            scsig = _sc(sig)
            scphi = _sc(tphi)
            if self.f > 0:
                # The general expression for tchi is
                #   tphi * scsig - sig * scphi
                # This involves cancellation if f > 0, so change to
                #   (tphi - sig) * (tphi + sig) / (tphi * scsig + sig * scphi)
                # To control overflow, write as (sigtphi = sig / tphi)
                #   (tphi - sig) * (1 + sigtphi) / (scsig + sigtphi * scphi)
                sigtphi = sig / tphi
                if sig < (tphi * _0_5):
                    t = tphi - sig
                else:
                    def _asinh_2(x):
                        return asinh(x) * _0_5
                    # Still possibly dangerous cancellation in tphi - sig.
                    # Write tphi - sig = (1 - e) * Dg(1, e)
                    #   Dg(x, y) = (g(x) - g(y)) / (x - y)
                    #   g(x) = sinh(x * atanh(sphi * x))
                    # Note sinh(atanh(sphi)) = tphi
                    # Turn the crank on divided differences, substitute
                    #   sphi = tphi / _sc(tphi)
                    #   atanh(x) = asinh(x / sqrt(1 - x^2))
                    e   = self._e
                    em1 = self._e2m1 / (_1_0 + e)
                    # assert em1 != 0
                    scb = self._scbeta(tphi)
                    scphib =  scphi / scb  # sec(phi) / sec(beta)
                    atphib = _asinh_2(tphi * e / scb)  # atanh(e * sphi)
                    atphi  = _asinh_2(tphi)  # atanh(sphi)
                    t      = _asinh_2(em1 * (tphi * scphib)) / em1
                    try:
                        Dg = Fsum(atphi, atphib, t, e * t)
                    except ValueError:  # Fsum(NAN) exception
                        Dg = _sum((atphi, atphib, t, e * t))
                    e *= atphib
                    t  = atphi - e
                    if t:  # sinh(0) == 0
                        Dg *= sinh(t) / t * cosh(atphi + e) * em1
                        t   = float(Dg)  # tphi - sig
                tchi = _over(t * (_1_0 + sigtphi),
                         scsig + scphi * sigtphi)  # if t else _0_0
            else:
                tchi =  tphi * scsig - sig * scphi

        n = _xkwds_get(diff_name, name=Phi.name)
        Chi = AuxChi(tchi, name=n).copyquadrant(Phi)

        if _xkwds_get(diff_name, diff=False):
            if isinf(tphi):  # PYCHOK np cover
                d = self._conformal_diff
            else:
                d = self.Parametric(Phi)._x_normalized
                if d:
                    d = _over(d, Chi._x_normalized) * \
                        _over(d, Phi._x_normalized) * self._e2m1
            Chi._diff = d
        # assrt Chi._AUX == Aux.CHI
        return Chi

    @Property_RO
    def _conformal_diff(self):  # PYCHOK no cover
        '''(INTERNAL) Constant I{Conformal} diff.
        '''
        e = self._e
        if self.f > 0:
            ss = sinh(asinh(self._e1) * e)
            d = _over(_1_0, _sc(ss) + ss)
        elif e:
            ss = sinh(-atan(e) * e)
            d = _sc(ss) - ss
        else:
            d = _1_0
        return d

    def convert(self, auxout, Zeta_d, exact=False):
        # Convert I{Auxiliary} or I{scalar} latitude
        Z = d = Zeta_d
        if isinstance(Z, AuxAngle):
            A, auxin = _AuxClass(auxout), Z._AUX
            if auxin == auxout:
                pass
            elif not (isfinite(Z.tan) and Z.tan):  # XXX
                Z = A(*Z._yx, aux=auxout, name=Z.name)
            elif exact:
                p = Aux.power(auxout, auxin)
                if p is None:
                    P = self._fromAux(Z)  # Phi
                    Z = self._toAux(auxout, P)
                    Z._iter = P.iteration
                else:
                    y, x = Z._yx
                    if p:
                        y *= pow(self._fm1, p)
                    Z = A(y, x, aux=auxout, name=Z.name)
            else:
                cs =  self._coeffs(auxout, auxin)
                yx =  Z._yx_normalized
                Z  =  A(*yx, aux=auxout, name=Z.name)
                # assert Z._yx == yx
                r  = _Clenshaw(True, Z, cs, self.ALorder)
                Z +=  AuxAngle.fromRadians(r, aux=auxout)
            # assert Z._AUX == auxout
            return Z

        elif isinstance(d, Degrees) or isscalar(d):
            Z  = AuxPhi.fromDegrees(d)
            d  = round((d - Z.toDegrees) / _360_0) * _360_0
            d += self.convert(auxout, Z, exact).toDegrees
            return Degrees(d, name=Aux.Greek(auxout))

        raise AuxError(auxout=auxout, Zeta_d=Zeta_d, exact=exact)

    def _Dq(self, tphi):
        # I{Divided Difference} of (q(1) - q(sphi)) / (1 - sphi).
        sphi = _sn(tphi)
        if tphi > 0:
            scphi = _sc(tphi)
            d = _1_0 / (scphi**2 * (_1_0 + sphi))  # XXX - sphi
            if d:
                # General expression for _Dq(1, sphi) is
                # atanh(e * d / (1 - e2 * sphi)) / (e * d) +
                # (1 + e2 * sphi) / ((1 - e2 * sphi * sphi) * e2m1)
                # with atanh(e * d / (1 - e2 * sphi)) =
                #      atanh(e * d * scphi / (scphi - e2 * tphi))
                e2m1, ed = self._e2m1, (self._e * d)
                if e2m1 and ed:
                    e2 = self._e2
                    if e2 > 0:  # assert self.f > 0
                        scb =  self._scbeta(tphi)
                        q   =  scphib = scphi / scb
                        q  *= (scphi + tphi * e2) / (e2m1 * scb)
                        q  +=  asinh(self._e1 * d * scphib) / ed
                    else:
                        s2  =  sphi * e2
                        q   = (_1_0 + s2) / ((_1_0 - sphi * s2) * e2m1)
                        q  += (atan2(ed, _1_0 - s2) / ed) if e2 < 0 else _1_0
                else:  # PYCHOK no cover
                    q = INF
            else:  # PYCHOK no cover
                q = self._2_e2m12
        else:  # not reached, open-coded in .Authalic
            q = _over(self._q - self._qf(tphi), _1_0 - sphi)
        return q

    @Property_RO
    def _e(self):  # unsigned, (1st) eccentricity
        return self.ellipsoid.e  # sqrt(fabs(self._e2))

    @Property_RO
    def _e1(self):
        return sqrt(fabs(self._e12))

    @Property_RO
    def _e12(self):
        return _over(self._e2, _1_0 - self._e2)

    @Property_RO
    def _e12p1(self):
        return _1_0 / self._e2m1

    @Property_RO
    def _e2(self):  # signed, (1st) eccentricity squared
        return self.ellipsoid.e2

    @Property_RO
    def _e2m1(self):  # 1 less 1st eccentricity squared
        return self.ellipsoid.e21  # == ._fm1**2

    @Property_RO
    def _e2m1_sq2(self):
        return self._e2m1 * sqrt(self._q * _0_5)

    @Property_RO
    def _2_e2m12(self):
        return _2_0 / self._e2m1**2

    @Property_RO
    def _Ef_fRG_a2b2_PI_4(self):
        E = self.ellipsoid
        return _Ef.fRG(E.a2, E.b2) / PI_4

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (L{Ellipsoid}).
        '''
        return self._E

    @Property_RO
    def f(self):
        '''Get the C{ellipsoid}'s flattening (C{scalar}).
        '''
        return self.ellipsoid.f

    flattening = f

    @Property_RO
    def _fm1(self):  # 1 - flattening
        return self.ellipsoid.f1

    def _fromAux(self, Zeta, **name):
        '''Convert I{Auxiliary} to I{Geographic} latitude.

           @arg Zeta: Auxiliary latitude (L{AuxAngle}).

           @return: Geographic latitude, I{Phi} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Zeta=Zeta)
        aux =  Zeta._AUX
        n   = _xkwds_get(name, name=Zeta.name)
        f   =  self._fromAuxCase.get(aux, None)
        if f is None:
            Phi = AuxPhi(NAN, name=n)
        elif callable(f):
            t   =  self._fm1
            t  *=  f(t)
            Phi = _Newton(t, Zeta, self._toZeta(aux), name=n)
        else:  # assert isscalar(f)
            y, x = Zeta._yx
            Phi  = AuxPhi(y / f, x, name=n)
        # assert Phi._AUX == Aux.PHI
        return Phi

    @Property_RO
    def _fromAuxCase(self):
        '''(INTERNAL) switch(auxin): ...
        '''
        return {Aux.AUTHALIC:         cbrt,
                Aux.CONFORMAL: _MODS.utily._passarg,
                Aux.GEOCENTRIC: self._e2m1,
                Aux.GEOGRAPHIC:      _1_0,
                Aux.PARAMETRIC: self._fm1,
                Aux.RECTIFYING:       sqrt}

    def Geocentric(self, Phi, **diff_name):
        '''Convert I{Geographic} to I{Geocentric} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Geocentric latitude, C{Phi} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        Theta = AuxTheta(Phi.y * self._e2m1, Phi.x,
                         name=_xkwds_get(diff_name, name=Phi.name))
        if _xkwds_get(diff_name, diff=False):
            Theta._diff = self._e2m1
        return Theta

    def Geodetic(self, Phi, **diff_name):  # PYCHOK no cover
        '''Convert I{Geographic} to I{Geodetic} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Geodetic latitude, C{Phi} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        return AuxPhi(Phi, name=_xkwds_get(diff_name, name=Phi.name))

    @Property_RO
    def _n(self):  # 3rd flattening
        return self.ellipsoid.n

    @Property_RO
    def _n2(self):
        return self._n**2

    def Parametric(self, Phi, **diff_name):
        '''Convert I{Geographic} to I{Parametric} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Parametric latitude, C{Beta} (L{AuxAngle}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        Beta = AuxBeta(Phi.y * self._fm1, Phi.x,
                       name=_xkwds_get(diff_name, name=Phi.name))
        if _xkwds_get(diff_name, diff=False):
            Beta._diff = self._fm1
        return Beta

    Reduced = Parametric

    @Property_RO
    def _q(self):  # constant _q
        q, f = self._e12p1, self.f
        if f:
            e  =  self._e
            q += (asinh(self._e1) if f > 0 else atan(e)) / e
        else:
            q += _1_0
        return q

    def _qf(self, tphi):
        # function _q: atanh(e * sphi) / e + sphi / (1 - (e * sphi)^2)
        scb  = self._scbeta(tphi)
        return self._atanhee(tphi) + (tphi / scb) * (_sc(tphi) / scb)

    def _qIntegrand(self, beta):
        # pbeta(beta) = integrate(q(beta), beta), with beta in radians
        #   q(beta) = (1-f) * (sin(xi) - sin(chi)) / cos(phi)
        #           = (1-f) * (cos(chi) - cos(xi)) / cos(phi) *
        #                     (cos(xi) + cos(chi)) / (sin(xi) + sin(chi))
        # Fit q(beta)/cos(beta) with Fourier transform
        #   q(beta)/cos(beta) = sum(c[k] * sin((2*k+1)*beta), k, 0, K-1)
        # then the integral is
        #   pbeta = sum(d[k] * cos((2*k+2)*beta), k, 0, K-1)
        # where
        #   d[k] = -1/(4*(k+1)) * (c[k] + c[k+1]) for k in 0..K-2
        #   d[K-1] = -1/(4*K) * c[K-1]
        Beta = AuxBeta.fromRadians(beta)
        if Beta.x:  # and self._fm1:
            Ax, _cv    = Aux, self.convert
            Phi        = _cv(Ax.PHI, Beta, exact=True)
            schi, cchi = _cv(Ax.CHI, Phi,  exact=True)._yx_normalized
            sxi,  cxi  = _cv(Ax.XI,  Phi,  exact=True)._yx_normalized
            r  = (sxi - schi) if fabs(schi) < fabs(cchi) else \
                 _over(_2cos2x(cchi, cxi), (sxi + schi) * _2_0)
            r *= _over(self._fm1, Phi._x_normalized * Beta._x_normalized)
        else:  # beta == PI_2, PI3_2, ...
            r  = _0_0  # XXX 0 avoids NAN summation exceptions
        return r

    def Rectifying(self, Phi, **diff_name):
        '''Convert I{Geographic} to I{Rectifying} latitude.

           @arg Phi: Geographic latitude (L{AuxAngle}).

           @return: Rectifying latitude, C{Mu} (L{AuxAngle}).
        '''
        Beta = self.Parametric(Phi)
        # assert Beta._AUX == Aux.BETA
        sb, cb = map(fabs, Beta._yx_normalized)
        a, ka, ka1 =      _1_0,  self._e2,  self._e2m1
        b, kb, kb1 = self._fm1, -self._e12, self._e12p1
        if self.f < 0:
            a,   b   = b,   a
            ka,  kb  = kb,  ka
            ka1, kb1 = kb1, ka1
            sb,  cb  = cb,  sb
        # now a, b = larger/smaller semiaxis
        # Beta measured from larger semiaxis
        # kb, ka = modulus-squared for distance from Beta = 0, pi/2
        # NB kb <= 0; 0 <= ka <= 1
        # sa = b*E(Beta, sqrt(kb))
        # sb = a*E(Beta',sqrt(ka))
        # 1 - ka * (1 - sb2) = 1 - ka + ka*sb2
        sb2 =  sb**2
        cb2 =  cb**2
        da2 =  ka1 + ka * sb2
        db2 = _1_0 - kb * sb2
        # DLMF Eq. 19.25.9
        my = b * sb * _Ef._RFRD(cb2, db2, _1_0, kb * sb2)
        # DLMF Eq. 19.25.10 with complementary angles
        mx = a * cb * (_Ef.fRF(sb2,  da2, _1_0) * ka1 +
            ka * cb2 * _Ef.fRD(sb2, _1_0,  da2, _3_0) * ka1 +
            ka * sb / sqrt(da2))
        # my + mx = 2*_Ef.fRG(a*a, b*b) = a*E(e) = b*E(i*e')
        # mr = a*E(e)*(2/pi) = b*E(i*e')*(2/pi)
        if self.f < 0:
            a,  b  = b,  a
            my, mx = mx, my
        mr = (my + mx) / PI_2
        if mr:
            my = sin(my / mr)
            mx = sin(mx / mr)  # XXX zero?
        else:  # zero Mu
            my, mx = _0_0, _1_0
        n  = _xkwds_get(diff_name, name=Phi.name)
        Mu =  AuxMu(my, mx,  # normalized
                    name=n).copyquadrant(Phi)

        if _xkwds_get(diff_name, diff=False):
            d, x = _0_0, Beta._x_normalized
            if x and mr:
                if Mu.x and Phi.x and not isinf(Phi.tan):
                    d = b / mr * (x / Mu.x)**2 \
                               * (x / Phi._x_normalized)
                else:
                    d = mr / a
            Mu._diff = self._fm1 * d
        return Mu

    def RectifyingRadius(self, exact=False):
        '''Get the I{Rectifying} radius.

           @arg exact: If C{True}, use the exact expression,
                       otherwise the I{Taylor} series.

           @return: Rectifying radius (L{Meter}, same units
                    as the ellipsoid axes).
        '''
        r = self._Ef_fRG_a2b2_PI_4 if exact else self._RectifyingR
        return Meter(r, name=self.RectifyingRadius.__name__)

    @Property_RO
    def _RectifyingR(self):
        m =  self.ALorder
        d = _MODS.karney._polynomial(self._n2, _RRCoeffs[m], 0, m // 2)
        return d * (self.a + self.b) * _0_5

    def _scbeta(self, tphi):
        return _sc(self._fm1 * tphi)

    def _toAux(self, auxout, Phi, **diff_name):
        '''Convert I{Geographic} to I{Auxiliary} latitude.

           @arg auxout: I{Auxiliary} kind (C{Aux.KIND}).
           @arg Phi: Geographic latitude (L{AuxLat}).

           @return: Auxiliary latitude, I{Eta} (L{AuxLat}).
        '''
        _xinstanceof(AuxAngle, Phi=Phi)
        # assert Phi._AUX == Aux.PHI
        n = _xkwds_get(diff_name, name=Phi.name)
        m = _toAuxCase.get(auxout, None)
        if m:  # callable
            A = m(self, Phi, **_xkwds(diff_name, name=n))
        elif auxout == Aux.GEODETIC:  # == GEOGRAPHIC
            A = AuxPhi(Phi, name=n)
        else:  # auxout?
            A = AuxPhi(NAN, name=n)
        # assert A._AUX == auxout
        return A

    def _toZeta(self, zetaux):
        '''Return a (lean) function to create C{AuxPhi(tphi)} and
           convert that into C{AuxAngle} of (fixed) kind C{zetaux}
           for use only inside the C{_Newton} loop.
        '''
        class _AuxPhy(AuxPhi):
            # lean C{AuxPhi} instance.
            # _diff = _1_0
            # _x    = _1_0

            def __init__(self, tphi):  # PYCHOK signature
                self._y = tphi

        m = _toAuxCase.get(zetaux, None)
        if m:  # callable

            def _toZeta(tphi):
                return m(self, _AuxPhy(tphi), diff=True)

        elif zetaux == Aux.GEODETIC:  # GEOGRAPHIC
            _toZeta = _AuxPhy

        else:  # zetaux?

            def _toZeta(unused):  # PYCHOK expected
                return _AuxPhy(NAN)

        return _toZeta


# switch(auxout): ...
_toAuxCase = {Aux.AUTHALIC:   AuxLat.Authalic,
              Aux.CONFORMAL:  AuxLat.Conformal,
              Aux.GEOCENTRIC: AuxLat.Geocentric,
              Aux.PARAMETRIC: AuxLat.Parametric,
              Aux.RECTIFYING: AuxLat.Rectifying}


def _Clenshaw(sinp, Zeta, cs, K):
    sz, cz = Zeta._yx  # isnormal
    # Evaluate sum(c[k] * sin((2*k+2) * Zeta)) if sinp else
    #          sum(c[k] * cos((2*k+2) * Zeta))
    x = _2cos2x(cz, sz)  # 2 * cos(2*Zeta)
    if isfinite(x):
        U0, U1 = Fsum(), Fsum()
        # assert len(cs) == K
        for r in _reverange(K):
            U1 -= U0 * x + cs[r]
            U1, U0 = U0, -U1
        # u0*f0(Zeta) - u1*fm1(Zeta)
        # f0  = sin(2*Zeta) if sinp else cos(2*Zeta)
        # fm1 = 0 if sinp else 1
        if sinp:
            U0 *= sz * cz * _2_0
        else:
            U0 *= x * _0_5
            U0 -= U1
        x = float(U0)
    return x


def _CXcoeffs(aL):  # PYCHOK in .auxilats.__main__
    '''(INTERNAL) Get the C{CX_4}, C{_6} or C{_8} coefficients.
    '''
    try:  # from pygeodesy.auxilats._CX_x import _coeffs_x as _coeffs
        _CX_x   = _DOT_(_MODS.auxilats.__name__, _UNDER_('_CX', aL))
        _coeffs = _MODS.getattr(_CX_x, _UNDER_('_coeffs', aL))
    except (AttributeError, ImportError, KeyError, TypeError) as x:
        raise AuxError(ALorder=aL, cause=x)
    # assert _coeffs.ALorder == aL
    # assert _coeffs.n == Aux.len(aL)
    return _coeffs


def _Newton(tphi, Zeta, _toZeta, **name):
    # Newton's method fro AuxLat._fromAux
    try:
        _lg2  = _log2
        _abs  =  fabs
        tz    = _abs(Zeta.tan)
        tphi  =  tz / tphi  # **)
        ltz   = _lg2(tz)    # **)
        ltphi = _lg2(tphi)  # **)
        ltmin =  min(ltphi, MIN_EXP)
        ltmax =  max(ltphi, MAX_EXP)
#       auxin =  Zeta._AUX
        s, n, __2  =  0, 3, _0_5  # n = i + 2
        _TOL, _xp2 = _EPSqrt, _exp2
        for i in range(1, _TRIPS):  # up to 1 Ki!
            # _toAux(auxin, AuxPhi(tphi), diff=True)
            Z = _toZeta(tphi)
            # assert Z._AUX == auxin
            t, s_ = Z.tan, s
            if t > tz:
                ltmax, s = ltphi, +1
            elif t < tz:
                ltmin, s = ltphi, -1
            else:
                break
            # derivative dtan(Z)/dtan(Phi)
            # to dlog(tan(Z))/dlog(tan(Phi))
            d = (ltz - _lg2(t)) * t  # **)
            if d:
                d = d / (Z.diff * tphi)  # **)
                ltphi += d
                tphi = _xp2(ltphi)
            if _abs(d) < _TOL:
                i += 1
                # _toAux(auxin, AuxPhi(tphi), diff=True)
                Z = _toZeta(tphi)
                tphi -= _over(Z.tan - tz, Z.diff)
                break
            if (i > n and (s * s_) < 0) or not ltmin < ltphi < ltmax:
                s, n  =  0, (i + 2)
                ltphi = (ltmin + ltmax) * __2
                tphi  = _xp2(ltphi)
        else:
            i = _TRIPS
        Phi = AuxPhi(tphi, **name).copyquadrant(Zeta)
        Phi._iter = i
    except (ValueError, ZeroDivisionError):  # **) zero t, tphi, tz or Z.diff
        Phi = AuxPhi(Zeta, **name)  # diff as-as
        Phi._iter = 0
    # assert Phi._AUX == Aux.PHI
    return Phi


_f, _u =  float, _Ufloats()
_1__f3 = -1 / _f(3)  # XXX +1 / _f(3)
_AR2Coeffs = {4: _u(4  / _f(315),   4 / _f(105),  4 / _f(15),  _1__f3),
              6: _u(4  / _f(1287),  4 / _f(693),  4 / _f(315),  4 / _f(105),
                    4  / _f(15),   _1__f3),
              8: _u(4  / _f(3315),  4 / _f(2145), 4 / _f(1287), 4 / _f(693),
                    4  / _f(315),   4 / _f(105),  4 / _f(15),  _1__f3)}
_RRCoeffs  = {4: _u(1  / _f(64),   _0_25),
              6: _u(1  / _f(256),   1 / _f(64),  _0_25),
              8: _u(25 / _f(16384), 1 / _f(256),  1 / _f(64),  _0_25)}  # PYCHOK used!
del _f, _u, _Ufloats, _1__f3
# assert set(_AR2Coeffs.keys()) == set(_RRCoeffs.keys())

__all__ += _ALL_DOCS(AuxLat)

# **) MIT License
#
# Copyright (C) 2023-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
