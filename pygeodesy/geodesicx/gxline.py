
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ class U{GeodesicLineExact
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicLineExact.html>}.

Class L{GeodesicLineExact} follows the naming, methods and return
values from class C{GeodesicLine} from I{Karney}'s Python U{geographiclib
<https://GeographicLib.SourceForge.io/1.52/python/index.html>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2023)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

# A copy of comments from Karney's C{GeodesicLineExact.cpp}:
#
# This is a reformulation of the geodesic problem.  The
# notation is as follows:
# - at a general point (no suffix or 1 or 2 as suffix)
#   - phi = latitude
#   - beta = latitude on auxiliary sphere
#   - omega = longitude on auxiliary sphere
#   - lambda = longitude
#   - alpha = azimuth of great circle
#   - sigma = arc length along great circle
#   - s = distance
#   - tau = scaled distance (= sigma at multiples of PI/2)
# - at northwards equator crossing
#   - beta = phi = 0
#   - omega = lambda = 0
#   - alpha = alpha0
#   - sigma = s = 0
# - a 12 suffix means a difference, e.g., s12 = s2 - s1.
# - s and c prefixes mean sin and cos

# from pygeodesy.basics import _xinstanceof  # from .karney
from pygeodesy.constants import NAN, _EPSmin, _0_0, _1_0, _180_0, _2__PI
from pygeodesy.fsums import _COMMASPACE_, fsumf_, fsum1f_
from pygeodesy.geodesicx.gxbases import _cosSeries, _GeodesicBase, \
                                        _sincos12, _sin1cos2
# from pygeodesy.interns import _COMMASPACE_  # from .fsums
from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS
from pygeodesy.karney import _around, _atan2d, Caps, _copysign, GDict, \
                             _fix90, _K_2_0, _norm2, _norm180, \
                             _sincos2, _sincos2d, _xinstanceof
from pygeodesy.props import Property_RO, _update_all
# from pygeodesy.streprs import pairs  # _MODS
from pygeodesy.utily import atan2d as _atan2d_reverse, sincos2

from math import atan2, cos, degrees, fabs, floor, radians, sin

__all__ = ()
__version__ = '23.08.19'

_glXs = []  # instances of C{[_]GeodesicLineExact} to be updated
# underflow guard, we require _TINY * EPS > 0, _TINY + EPS == EPS
_TINY = _EPSmin
# assert (_TINY * EPS) > 0 and (_TINY + EPS) == EPS


def _update_glXs(gX):  # see GeodesicExact.C4order and -._ef_reset_k2
    '''(INTERNAL) Zap cached/memoized C{Property[_RO]}s of
       any L{GeodesicLineExact} instances tied to the given
       L{GeodesicExact} instance B{C{gX}}.
    '''
    _xinstanceof(_MODS.geodesicx.GeodesicExact, gX=gX)
    for glX in _glXs:  # PYCHOK use weakref?
        if glX._gX is gX:
            _update_all(glX)


class _GeodesicLineExact(_GeodesicBase):
    '''(INTERNAL) Base class for L{GeodesicLineExact}.
    '''
    _a13   = _s13 = NAN
    _azi1  = _0_0
    _cchi1 =  NAN
    _dn1   =  NAN
    _gX    =  None  # Exact only
    _k2    =  NAN
    _lat1  = _lon1 = _0_0
    _salp0 = _calp0 = NAN
    _salp1 = _calp1 = NAN
    _somg1 = _comg1 = NAN
    _ssig1 = _csig1 = NAN

    def __init__(self, gX, lat1, lon1, azi1, caps, _debug, *salp1_calp1, **name):
        '''(INTERNAL) New C{[_]GeodesicLineExact} instance.
        '''
        _xinstanceof(_MODS.geodesicx.GeodesicExact, gX=gX)
        Cs = Caps
        if _debug:  # PYCHOK no cover
            self._debug |= _debug & Cs._DEBUG_ALL
            # _CapsBase.debug._update(self)
        if salp1_calp1:
            salp1, calp1 = salp1_calp1
        else:
            azi1 = _norm180(azi1)
            # guard against salp0 underflow,
            # also -0 is converted to +0
            salp1, calp1 = _sincos2d(_around(azi1))
        if name:
            self.name = name

        self._gX    = gX  # GeodesicExact only
        self._lat1  = lat1 = _fix90(lat1)
        self._lon1  = lon1
        self._azi1  = azi1
        self._salp1 = salp1
        self._calp1 = calp1
        # allow lat, azimuth and unrolling of lon
        self._caps  = caps | Cs._LINE

        sbet1, cbet1 = gX._sinf1cos2d(_around(lat1))
        self._dn1 = gX._dn(sbet1, cbet1)
        # Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0), with alp0
        # in [0, pi/2 - |bet1|].  Alt: calp0 = hypot(sbet1, calp1 * cbet1),
        # but the following is slightly better, consider the case salp1 = 0.
        self._salp0, self._calp0 = _sin1cos2(salp1, calp1, sbet1, cbet1)
        self._k2 = self._calp0**2 * gX.ep2
        # Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
        # sig = 0 is nearest northward crossing of equator.
        # With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
        # With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
        # With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
        # Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
        # With alp0 in (0, pi/2], quadrants for sig and omg coincide.
        # No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
        # With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
        self._somg1 = sbet1 * self._salp0
        self._comg1 = c = (cbet1 * calp1) if (sbet1 or calp1) else _1_0
        # Without normalization we have schi1 = somg1.
        self._cchi1 = gX.f1 * self._dn1 * c
        self._ssig1, self._csig1 = _norm2(sbet1, c)  # sig1 in (-pi, pi]
        # _norm2(somg1, comg1)  # no need to normalize!
        # _norm2(schi1?, cchi1)  # no need to normalize!
        if not (caps & Cs.LINE_OFF):
            _glXs.append(self)
        # no need to pre-compute other attrs based on _Caps.X.  All are
        # Property_RO's, computed once and cached/memoized until reset
        # when C4order is changed or Elliptic function reset is invoked.

    def __del__(self):  # XXX use weakref?
        if _glXs:  # may be empty or None
            try:  # PYCHOK no cover
                _glXs.remove(self)
            except (TypeError, ValueError):
                pass
        self._gX = None
        # _update_all(self)  # throws TypeError during Python 2 cleanup

    def _update(self, updated, *attrs, **unused):
        if updated:
            _update_all(self, *attrs)

    @Property_RO
    def a1(self):
        '''Get the I{equatorial arc} (C{degrees}), the arc length between
           the northward equatorial crossing and the first point.
        '''
        return _atan2d(self._ssig1, self._csig1)  # or NAN

    equatorarc = a1

    @Property_RO
    def a13(self):
        '''Get the arc length to reference point 3 (C{degrees}).

           @see: Methods L{Arc} and L{SetArc}.
        '''
        return self._a13

    def Arc(self):
        '''Return the arc length to reference point 3 (C{degrees} or C{NAN}).

           @see: Method L{SetArc} and property L{a13}.
        '''
        return self.a13

    def ArcPosition(self, a12, outmask=Caps.STANDARD):
        '''Find the position on the line given B{C{a12}}.

           @arg a12: Spherical arc length from the first point to the
                     second point (C{degrees}).
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and arc length C{a12} always included,
                    except when C{a12=NAN}.

           @note: By default, C{B{outmask}=STANDARD}, meaning thc C{lat1},
                  C{lon1}, C{azi1}, C{lat2}, C{lon2}, C{azi2}, C{s12} and
                  C{a12} entries are returned, except when C{a12=NAN}.
        '''
        return self._GDictPosition(True, a12, outmask)

    @Property_RO
    def azi0(self):
        '''Get the I{equatorial azimuth}, the azimuth of this geodesic line
           as it crosses the equator in a northward direction (C{degrees90}).
        '''
        return _atan2d(*self.azi0_sincos2)  # or NAN

    equatorazimuth = azi0

    @Property_RO
    def azi0_sincos2(self):
        '''Get the sine and cosine of the I{equatorial azimuth} (2-tuple C{(sin, cos)}).
        '''
        return self._salp0, self._calp0

    @Property_RO
    def azi1(self):
        '''Get the azimuth at the first point (compass C{degrees}).
        '''
        return self._azi1

    @Property_RO
    def azi1_sincos2(self):
        '''Get the sine and cosine of the first point's azimuth (2-tuple C{(sin, cos)}).
        '''
        return self._salp1, self._calp1

    @Property_RO
    def _B41(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return _cosSeries(self._C4a, self._ssig1, self._csig1)

    @Property_RO
    def _C4a(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self.geodesic._C4f_k2(self._k2)

    @Property_RO
    def _caps_DISTANCE_IN(self):
        '''(INTERNAL) Get C{Caps.DISTANCE_IN} and C{_OUT}.
        '''
        return self.caps & Caps._DISTANCE_IN_OUT

    @Property_RO
    def _D0k2(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.cD * _2__PI * self._k2

    @Property_RO
    def _D1(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.deltaD(self._ssig1, self._csig1, self._dn1)

    def Distance(self):
        '''Return the distance to reference point 3 (C{meter} or C{NAN}).

           @see: Method L{SetDistance} and property L{s13}.
        '''
        return self.s13

    @Property_RO
    def _E0b(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.cE * _2__PI * self.geodesic.b

    @Property_RO
    def _E1(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.deltaE(self._ssig1, self._csig1, self._dn1)

    @Property_RO
    def _eF(self):
        '''(INTERNAL) Cached/memoized C{Elliptic} function.
        '''
        # see .gx.GeodesicExact._ef_reset_k2
        return _MODS.elliptic.Elliptic(k2=-self._k2, alpha2=-self.geodesic.ep2)

    def _GDictPosition(self, arcmode, s12_a12, outmask=Caps.STANDARD):  # MCCABE 17
        '''(INTERNAL) Generate a new position along the geodesic.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and arc length C{a12} always included,
                    except when C{a12=NAN}.
        '''

        r = GDict(a12=NAN, s12=NAN)  # note both a12 and s12, always
        if not (arcmode or self._caps_DISTANCE_IN):  # PYCHOK no cover
            return r  # Uninitialized or impossible distance requested

        Cs = Caps
        if self._debug:  # PYCHOK no cover
            outmask |= self._debug & Cs._DEBUG_DIRECT_LINE
        outmask &= self._caps & Cs._OUT_MASK

        eF = self._eF
        gX = self.geodesic  # ._gX

        if arcmode:
            # s12_a12 is spherical arc length
            E2 = _0_0
            sig12 = radians(s12_a12)
            if _K_2_0:
                ssig12, csig12 = sincos2(sig12)  # utily, no NEG0
            else:  # PYCHOK no cover
                a  = fabs(s12_a12)  # 0 <= fabs(_remainder(s12_a12, _180_0)) <= 90
                a -= floor(a / _180_0) * _180_0  # 0 <= 0 < 180
                ssig12 = _0_0 if a ==  0 else sin(sig12)
                csig12 = _0_0 if a == 90 else cos(sig12)
        else:  # s12_a12 is distance
            t = s12_a12 / self._E0b
            s, c = _sincos2(t)  # tau12
            # tau2 = tau1 + tau12
            E2 = -eF.deltaEinv(*_sincos12(-s, c, *self._stau1_ctau1))
            sig12 = fsum1f_(self._E1, -E2, t)  # == t - (E2 - E1)
            ssig12, csig12 = _sincos2(sig12)

        salp0, calp0 = self._salp0, self._calp0
        ssig1, csig1 = self._ssig1, self._csig1

        # sig2 = sig1 + sig12
        ssig2, csig2 = _sincos12(-ssig12, csig12, ssig1, csig1)
        dn2 = eF.fDelta(ssig2, csig2)
        # sin(bet2) = cos(alp0) * sin(sig2) and
        #    cbet2  = hypot(salp0, calp0 * csig2).  Alt:
        #    cbet2  = hypot(csig2, salp0 * ssig2)
        sbet2, cbet2 = _sin1cos2(calp0, salp0, csig2, ssig2)
        if cbet2 == 0:  # salp0 = 0, csig2 = 0, break degeneracy
            cbet2 = csig2 = _TINY
        # tan(alp0) = cos(sig2) * tan(alp2)
        salp2 = salp0
        calp2 = calp0 * csig2  # no need to normalize

        if (outmask & Cs.DISTANCE):
            if arcmode:  # or f_0_01
                E2 = eF.deltaE(ssig2, csig2, dn2)
                # AB1 = _E0 * (E2 - _E1)
                # s12 = _b * (_E0 * sig12 + AB1)
                #     = _b * _E0 * (sig12 + (E2 - _E1))
                #     = _b * _E0 * (E2 - _E1 + sig12)
                s12 = self._E0b * fsum1f_(E2, -self._E1, sig12)
            else:
                s12 = s12_a12
            r.set_(s12=s12)

        if (outmask & Cs._DEBUG_DIRECT_LINE):  # PYCHOK no cover
            r.set_(sig12=sig12, dn2=dn2, b=gX.b, e2=gX.e2, f1=gX.f1,
                   E0b=self._E0b, E1=self._E1, E2=E2, eFk2=eF.k2, eFa2=eF.alpha2)

        if (outmask & Cs.LONGITUDE):
            schi1 = self._somg1
            cchi1 = self._cchi1
            schi2 = ssig2 * salp0
            cchi2 = gX.f1 * dn2 * csig2  # schi2 = somg2 without normalization
            lam12 = salp0 * self._H0e2_f1 * fsum1f_(eF.deltaH(ssig2, csig2, dn2),
                                                    -self._H1, sig12)
            if (outmask & Cs.LONG_UNROLL):
                t = _copysign(_1_0, salp0)  # east-going?
                tchi1 = t * schi1
                tchi2 = t * schi2
                chi12 = t * fsum1f_(atan2(ssig1, csig1), -atan2(ssig2, csig2),
                                    atan2(tchi2, cchi2), -atan2(tchi1, cchi1), sig12)
                lon2 = self.lon1 + degrees(chi12 - lam12)
            else:
                chi12 = atan2(*_sincos12(schi1, cchi1, schi2, cchi2))
                lon2 = _norm180(self._lon1_norm180 + _norm180(degrees(chi12 - lam12)))
            r.set_(lon2=lon2)
            if (outmask & Cs._DEBUG_DIRECT_LINE):  # PYCHOK no cover
                r.set_(ssig2=ssig2, chi12=chi12, H0e2_f1=self._H0e2_f1,
                       csig2=csig2, lam12=lam12, H1=self._H1)

        if (outmask & Cs.LATITUDE):
            r.set_(lat2=_atan2d(sbet2, gX.f1 * cbet2))

        if (outmask & Cs.AZIMUTH):
            r.set_(azi2=_atan2d_reverse(salp2, calp2, reverse=outmask & Cs.REVERSE2))

        if (outmask & Cs._REDUCEDLENGTH_GEODESICSCALE):
            dn1 = self._dn1
            J12 = self._D0k2 * fsumf_(eF.deltaD(ssig2, csig2, dn2), -self._D1, sig12)
            if (outmask & Cs._DEBUG_DIRECT_LINE):  # PYCHOK no cover
                r.set_(ssig1=ssig1, dn1=dn1, D0k2=self._D0k2,
                       csig1=csig1, J12=J12, D1=self._D1)
            if (outmask & Cs.REDUCEDLENGTH):
                # Add parens around (csig1 * ssig2) and (ssig1 * csig2) to
                # ensure accurate cancellation in the case of coincident points.
                r.set_(m12=gX.b * fsum1f_(dn2 * (csig1 * ssig2),
                                         -dn1 * (ssig1 * csig2),
                                         -J12 * (csig1 * csig2)))
            if (outmask & Cs.GEODESICSCALE):
                t = self._k2 * (ssig2 - ssig1) * (ssig2 + ssig1) / (dn2 + dn1)
                r.set_(M12=csig12 + ssig1 * (t * ssig2 - csig2 * J12) / dn1,
                       M21=csig12 - ssig2 * (t * ssig1 - csig1 * J12) / dn2)

        if (outmask & Cs.AREA):
            A4 = salp0 * calp0
            if A4:
                # tan(alp) = tan(alp0) * sec(sig)
                # tan(alp2-alp1) = (tan(alp2) - tan(alp1)) / (tan(alp2) * tan(alp1) + 1)
                # = calp0 * salp0 * (csig1 - csig2) / (salp0^2 + calp0^2 * csig1 * csig2)
                # If csig12 > 0, write
                #   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
                # else
                #   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
                # No need to normalize
                salp12 = (((ssig12 * csig1 / (_1_0 + csig12) + ssig1) * ssig12) if csig12 > 0 else
                                    (csig1 * (_1_0 - csig12) + ssig1  * ssig12)) * A4
                calp12 = salp0**2 + calp0**2 * csig1 * csig2
                A4 *=  gX._e2a2
                B41 =  self._B41
                B42 = _cosSeries(self._C4a, ssig2, csig2)
                S12 = (B42 - B41) * A4
            else:
                S12 = A4 = B41 = B42 = _0_0
                # alp12 = alp2 - alp1, used in atan2 so no need to normalize
                salp12, calp12 = _sincos12(self._salp1, self._calp1, salp2, calp2)
                # We used to include some patch up code that purported to deal
                # with nearly meridional geodesics properly.  However, this turned
                # out to be wrong once salp1 = -0 was allowed (via InverseLine).
                # In fact, the calculation of {s,c}alp12 was already correct
                # (following the IEEE rules for handling signed zeros).  So,
                # the patch up code was unnecessary (as well as dangerous).
            if (outmask & Cs._DEBUG_DIRECT_LINE):  # PYCHOK no cover
                r.set_(salp12=salp12, salp0=salp0, B41=B41, A4=A4,
                       calp12=calp12, calp0=calp0, B42=B42, c2=gX.c2)
            S12 += gX.c2 * atan2(salp12, calp12)
            r.set_(S12=S12)

        r.set_(a12=s12_a12 if arcmode else degrees(sig12),
               lat1=self.lat1,  # == _fix90(lat1)
               lon1=self.lon1 if (outmask & Cs.LONG_UNROLL) else self._lon1_norm180,
               azi1=_norm180(self.azi1))
        return r

    def _GenPosition(self, arcmode, s12_a12, outmask):
        '''(INTERNAL) Generate a new position along the geodesic.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2,
                                      s12, m12,  M12,  M21, S12)}.
        '''
        r = self._GDictPosition(arcmode, s12_a12, outmask)
        return r.toDirect9Tuple()

    def _GenSet(self, arcmode, s13_a13):
        '''(INTERNAL) Aka C++ C{GenSetDistance}.
        '''
        if arcmode:
            self.SetArc(s13_a13)
        else:
            self.SetDistance(s13_a13)
        return self  # for gx.GeodesicExact.InverseLine and -._GenDirectLine

    @Property_RO
    def geodesic(self):
        '''Get the I{exact} geodesic (L{GeodesicExact}).
        '''
        assert isinstance(self._gX, _MODS.geodesicx.GeodesicExact)
        return self._gX

    @Property_RO
    def _H0e2_f1(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.cH * _2__PI * self.geodesic._e2_f1

    @Property_RO
    def _H1(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self._eF.deltaH(self._ssig1, self._csig1, self._dn1)

    @Property_RO
    def lat1(self):
        '''Get the latitude of the first point (C{degrees}).
        '''
        return self._lat1

    @Property_RO
    def lon1(self):
        '''Get the longitude of the first point (C{degrees}).
        '''
        return self._lon1

    @Property_RO
    def _lon1_norm180(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return _norm180(self._lon1)

    def Position(self, s12, outmask=Caps.STANDARD):
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second (C{meter}).
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and arc length C{a12} always included,
                    except when C{a12=NAN}.

           @note: By default, C{B{outmask}=STANDARD}, meaning thc C{lat1},
                  C{lon1}, C{azi1}, C{lat2}, C{lon2}, C{azi2}, C{s12} and
                  C{a12} entries are returned, except when C{a12=NAN}.

           @note: This L{GeodesicLineExact} instance must have been
                  constructed with capability C{Caps.DISTANCE_IN} set.
        '''
        return self._GDictPosition(False, s12, outmask)

    @Property_RO
    def s13(self):
        '''Get the distance to reference point 3 (C{meter} or C{NAN}).

           @see: Methods L{Distance} and L{SetDistance}.
        '''
        return self._s13

    def SetArc(self, a13):
        '''Set reference point 3 in terms relative to the first point.

           @arg a13: Spherical arc length from the first to the reference
                     point (C{degrees}).

           @return: The distance C{s13} (C{meter}) between the first and
                    the reference point or C{NAN}.
        '''
        self._a13 = a13
        self._s13 = s13 = self._GDictPosition(True, a13, Caps.DISTANCE).s12
        _update_all(self)
        return s13

    def SetDistance(self, s13):
        '''Set reference point 3 in terms relative to the first point.

           @arg s13: Distance from the first to the reference point (C{meter}).

           @return: The arc length C{a13} (C{degrees}) between the first
                    and the reference point or C{NAN}.
        '''
        self._s13 = s13
        self._a13 = a13 = self._GDictPosition(False, s13, 0).a12
        _update_all(self)
        return a13  # NAN for GeodesicLineExact without Cap.DISTANCE_IN

    @Property_RO
    def _stau1_ctau1(self):
        '''(INTERNAL) Cached/memoized.
        '''
        s, c = _sincos2(self._E1)
        # tau1 = sig1 + B11
        return _sincos12(-s, c, self._ssig1, self._csig1)
        # unnecessary because Einv inverts E
        # return -self._eF.deltaEinv(stau1, ctau1)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{GeodesicLineExact} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: C{GeodesicLineExact} (C{str}).
        '''
        d = dict(geodesic=self.geodesic,
                 lat1=self.lat1, lon1=self.lon1, azi1=self.azi1,
                 a13=self.a13, s13=self.s13)
        return sep.join(_MODS.streprs.pairs(d, prec=prec))


__all__ += _ALL_DOCS(_GeodesicLineExact)

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
