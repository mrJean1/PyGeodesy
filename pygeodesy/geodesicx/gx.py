
# -*- coding: utf-8 -*-

u'''A Python version of I{Karney}'s C++ class U{GeodesicExact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}.

Classes L{GeodesicExact} and L{GeodesicLineExact} follow the naming,
methods and return values from I{Karney}s' Python classes C{Geodesic}
and C{GeodesicLine}, respectively.

Copyright (C) Charles Karney (2012-2021) <Charles@Karney.com>
and licensed under the MIT/X11 License.  For more information,
see U{GeographicLib<https://GeographicLib.SourceForge.io>}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

# A copy of comments from Karney's C{GeodesicExact.cpp}:
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

from pygeodesy.basics import copysign0, _xinstanceof, _xor, unsigned0
from pygeodesy.fmath import cbrt, fsum_, hypot, norm2
from pygeodesy.geodesicx.gxbases import _ALL_DOCS, Caps, _coSeries, \
                                        _GeodesicBase, _polynomial, \
                                        _sincos12, _TINY, _xnC4
from pygeodesy.geodesicx.gxline import _GeodesicLineExact, pairs, \
                                       _update_glXs
from pygeodesy.interns import EPS, EPS0, EPS02, MANT_DIG, NAN, NN, PI, PI_2, \
                             _COMMASPACE_, _convergence_, _EPSqrt, _no_, \
                             _0_0, _0_001, _0_01, _0_1, _0_5, _1_0, \
                             _N_1_0, _2_0, _3_0, _4_0, _6_0, _8_0, _16_0, \
                             _90_0, _180_0, _1000_0
from pygeodesy.karney import _around, _ellipsoid, _EWGS84, GDict, GeodesicError, \
                             _diff182, _fix90, _norm180, Property, Property_RO
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
# from pygeodesy.props import Property, Property_RO  # from .karney
# from pygeodesy.streprs import pairs  # from .geodesicx.gxline
from pygeodesy.utily import atan2d, sincos2, sincos2d, unroll180, wrap360

from math import atan2, cos, degrees, radians, sqrt

__all__ = ()
__version__ = '21.11.30'

_MAXIT1  = 20
_MAXIT2  = 10 + _MAXIT1 + MANT_DIG  # MANT_DIG == C++ digits
_SQRT2_2 = -0.7071  # negative sqrt(2) / 2
_1_75    =  1.75

# increased multiplier in defn of _TOL1 from 100 to 200 to fix Inverse
# case 52.784459512564 0 -52.784459512563990912 179.634407464943777557
# which otherwise failed for Visual Studio 10 (Release and Debug)
_TOL0 =  EPS
_TOL1 = _TOL0 * -200  # negative
_TOL2 = _EPSqrt  # == sqrt(_TOL0)
_TOL3 = _TOL2 * _0_1
_TOLb = _TOL2 * _TOL0  # Check on bisection interval
_THR1 = _TOL2 * _1000_0 + _1_0

_TINY3  = _TINY *  _3_0
_TOL08  = _TOL0 *  _8_0
_TOL016 = _TOL0 * _16_0


def _eTOL2(f):
    # Using the auxiliary sphere solution with dnm computed at
    # (bet1 + bet2) / 2, the relative error in the azimuth
    # consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
    # (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.

    # For a given f and sig12, the max error occurs for lines
    # near the pole.  If the old rule for computing dnm = (dn1
    # + dn2)/2 is used, then the error increases by a factor of
    # 2.)  Setting this equal to epsilon gives sig12 = etol2.

    # Here 0.1 is a safety factor (error decreased by 100) and
    # max(0.001, abs(f)) stops etol2 getting too large in the
    # nearly spherical case.
    t = min(_1_0, _1_0 - f * _0_5) * max(_0_001, abs(f)) * _0_5
    return _TOL3 / (sqrt(t) if t > EPS02 else EPS0)


class _PDict(GDict):
    '''(INTERNAL) Parameters passed around in C{._GDictInverse} and
       optionally returned when C{GeodesicExact.debug} is C{True}.
    '''
    def setsigs(self, ssig1, csig1, ssig2, csig2):
        '''Update the C{sig1} and C{sig2} parameters.
        '''
        self.set_(ssig1=ssig1, csig1=csig1, sncndn1=(ssig1, csig1, self.dn1),  # PYCHOK dn1
                  ssig2=ssig2, csig2=csig2, sncndn2=(ssig2, csig2, self.dn2))  # PYCHOK dn2

    def toGDict(self):
        '''Return as C{GDict} without attrs C{sncndn1} and C{sncndn2}.
        '''
        def _rest(sncndn1=None, sncndn2=None, **rest):  # PYCHOK unused
            return GDict(rest)

        return _rest(**self)


class GeodesicExact(_GeodesicBase):
    '''A pure Python version of I{Karney}'s C++ class U{GeodesicExact
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>},
       modeled after I{Karney}'s Python class U{Geodesic<https://GeographicLib.SourceForge.io/
       html/python/code.html#module-geographiclib.geodesic>}.
    '''
    _E   = _EWGS84
    _nC4 =  30  # default C4Order

    def __init__(self, a_ellipsoid=_EWGS84, f=None, name=NN, C4Order=None):
        '''New L{GeodesicExact} instance.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum}) or
                             the equatorial radius of the ellipsoid (C{scalar},
                             conventionally in C{meter}), see B{C{f}}.
           @arg f: The flattening of the ellipsoid (C{scalar}) if B{C{a_ellipsoid}}
                   is specified as C{scalar}.
           @kwarg name: Optional name (C{str}).
           @kwarg C4Order: Optional series expansion order (C{int}), see property
                           L{C4Order}, default C{30}.

           @raise GeodesicError: Invalid B{C{C4Order}}.
        '''
        if a_ellipsoid not in (GeodesicExact._E, None):
            self._E = _ellipsoid(a_ellipsoid, f, name=name)

        if name:
            self.name = name
        if C4Order:  # XXX private copy, always?
            self.C4Order = C4Order

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    def ArcDirect(self, lat1, lon1, azi1, a12, outmask=Caps.STANDARD):
        '''Solve the I{Direct} geodesic problem in terms of (spherical) arc length.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @arg a12: Arc length between the points (C{degrees}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and arc length C{a12} always included.

           @see: C++ U{GeodesicExact.ArcDirect
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.ArcDirect<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return self._GDictDirect(lat1, lon1, azi1, True, a12, outmask)

    def ArcDirectLine(self, lat1, lon1, azi1, a12, caps):
        '''Define a L{GeodesicLineExact} in terms of the I{direct} geodesic problem and as arc length.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @arg a12: Arc length between the points (C{degrees}), can be negative.
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineExact} instance
                        should possess, i.e., which quantities can be
                        returned by calls to L{GeodesicLineExact.Position}
                        and L{GeodesicLineExact.ArcPosition}.

           @return: A L{GeodesicLineExact} instance.

           @note: The third point of the L{GeodesicLineExact} is set to correspond
                  to the second point of the I{Inverse} geodesic problem.

           @note: Latitude B{C{lat1}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.ArcDirectLine
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.ArcDirectLine<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return self._GenDirectLine(lat1, lon1, azi1, True, a12, caps)

    def Area(self, polyline=False, name=NN):
        '''Set up a L{GeodesicAreaExact} to compute area and
           perimeter of a polygon.

           @kwarg polyline: If C{True} perimeter only, otherwise
                            area and perimeter (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: A L{GeodesicAreaExact} instance.

           @note: The B{C{debug}} setting is passed as C{verbose}
                  to the returned L{GeodesicAreaExact} instance.
        '''
        from pygeodesy.geodesicx.gxarea import GeodesicAreaExact
        gaX = GeodesicAreaExact(self, polyline=polyline,
                                      name=name or self.name)
        if self.debug:
            gaX.verbose = True
        return gaX

    Polygon = Area  # for C{geographiclib} compatibility

    @Property_RO
    def b(self):
        '''Get the ellipsoid's I{polar} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.b

    @Property_RO
    def c2x(self):
        '''Get the ellipsoid's I{authalic} earth radius I{squared} (C{meter**2}).
        '''
        # The Geodesic class substitutes atanh(sqrt(e2)) for asinh(sqrt(ep2))
        # in the definition of _c2.  The latter is more accurate for very
        # oblate ellipsoids (which the Geodesic class does not handle).  Of
        # course, the area calculation in GeodesicExact is still based on a
        # series and only holds for moderately oblate (or prolate) ellipsoids.
        return self.ellipsoid.c2x

    c2 = c2x  # in this particular case

    def C4f(self, eps):
        '''Evaluate the C{C4x} coefficients for B{C{eps}}.

           @arg eps: Polynomial factor (C{float}).

           @return: C{C4Order}-Tuple of C{C4x(B{eps})} coefficients.
        '''
        def _c4(nC4, C4x, _p_eps, eps):
            i, e = 0, _1_0
            for r in range(nC4, 0, -1):
                j = i + r
                yield _p_eps(eps, C4x, i, j) * e
                e *= eps
                i = j
            # assert i == (nC4 * (nC4 + 1)) // 2

        return tuple(_c4(self._nC4, self._C4x,
                        _polynomial, eps))

    def _C4f_k2(self, k2):  # in ._GDictInverse and gxline._GeodesicExactLine._C4a
        '''(INTERNAL) Compute C{eps} from B{C{k2}} and invoke C{C4f}.
        '''
        return self.C4f(k2 / fsum_(_2_0, sqrt(k2 + _1_0) * _2_0, k2))

    @Property
    def C4Order(self):
        '''Get the series expansion order (C{int}, 24, 27 or 30).
        '''
        return self._nC4

    @C4Order.setter  # PYCHOK .setter
    def C4Order(self, order):
        '''Set the series expansion order.

           @arg order: New order (C{int}, 24, 27 or 30).

           @raise GeodesicError: Invalid B{C{order}}.
        '''
        _xnC4(C4Order=order)
        if self._nC4 != order:
            GeodesicExact._C4x._update(self)
            _update_glXs(self)  # zap cached/memoized _GeodesicLineExact attrs
        self._nC4 = order

    @Property_RO
    def _C4x(self):
        '''Get this ellipsoid's C{C4} coefficients, I{cached} tuple.

           @see: Property L{C4Order}.
        '''
        # see C4coeff() in GeographicLib.src.GeodesicExactC4.cpp
        def _C4(nC4, coeffs, _p_n, n):
            i = 0
            for r in range(nC4 + 1, 1, -1):
                for j in range(1, r):
                    j = i + j  # (j - i - 1) order of polynomial
                    yield _p_n(n, coeffs, i, j) / coeffs[j]
                    i = j + 1
            # assert i == len(cs) == (nC4 * (nC4 + 1) * (nC4 + 5)) // 6

        return tuple(_C4(self._nC4, self._coeffs(self._nC4),
                        _polynomial, self.n))  # 3rd flattening

    def _coeffs(self, nC4):
        '''(INTERNAL) Get the series coefficients.
        '''
        if nC4 == 30:
            from pygeodesy.geodesicx._C4_30 import _coeffs_30 as _coeffs
        elif nC4 == 27:
            from pygeodesy.geodesicx._C4_27 import _coeffs_27 as _coeffs
        elif nC4 == 24:
            from pygeodesy.geodesicx._C4_24 import _coeffs_24 as _coeffs
        else:  # PYCHOK no cover
            _coeffs = ()
        n = _xnC4(nC4=nC4)
        if len(_coeffs) != n:  # double check
            raise GeodesicError(_coeffs=len(_coeffs), _xnC4=n, nC4=nC4)
        return _coeffs

    def Direct(self, lat1, lon1, azi1, s12, outmask=Caps.STANDARD):
        '''Solve the I{Direct} geodesic problem

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @arg s12: Distance between the points (C{meter}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and distance C{s12} always included.

           @see: C++ U{GeodesicExact.Direct
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Direct<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return self._GDictDirect(lat1, lon1, azi1, False, s12, outmask)

    def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
        '''Return the destination lat, lon and reverse azimuth
           (final bearing) in C{degrees}.

           @return: L{Destination3Tuple}C{(lat, lon, final)}.
        '''
        r = self._GDictDirect(lat1, lon1, azi1, False, s12, Caps._AZIMUTH_LATITUDE_LONGITUDE)
        return Destination3Tuple(r.lat2, r.lon2, r.azi2)

    def DirectLine(self, lat1, lon1, azi1, s12, caps=Caps.STANDARD):
        '''Define a L{GeodesicLineExact} in terms of the I{direct} geodesic problem and as distance.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @arg s12: Distance between the points (C{meter}), can be negative.
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineExact} instance
                        should possess, i.e., which quantities can be
                        returned by calls to L{GeodesicLineExact.Position}.

           @return: A L{GeodesicLineExact} instance.

           @note: The third point of the L{GeodesicLineExact} is set to correspond
                  to the second point of the I{Inverse} geodesic problem.

           @note: Latitude B{C{lat1}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.DirectLine
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.DirectLine<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return self._GenDirectLine(lat1, lon1, azi1, False, s12, caps)

    def _dn(self, sbet, cbet):  # in gxline._GeodesicLineExact.__init__
        '''(INTERNAL) Helper.
        '''
        return sqrt(_1_0 + self.ep2 * sbet**2) if self.f >= 0 else (
               sqrt(_1_0 - self.e2  * cbet**2) /  self.f1)

    @Property_RO
    def e2(self):
        '''Get the ellipsoid's I{(1st) eccentricity squared} (C{float}), M{f * (2 - f)}.
        '''
        return self.ellipsoid.e2

    @Property_RO
    def _e2a2(self):
        '''(INTERNAL) Cache M{E.e2 * E.a2}.
        '''
        return self.e2 * self.ellipsoid.a2

    @Property_RO
    def _e2_f1(self):
        '''(INTERNAL) Cache M{E.e2 * E.f1}.
        '''
        return self.e2 / self.f1

    @Property_RO
    def _eF(self):
        '''(INTERNAL) Get the elliptic function, aka C{.E}.
        '''
        from pygeodesy.elliptic import Elliptic
        return Elliptic(k2=-self.ep2)

    def _eF_reset_e2_f1_cH(self, x, y):
        '''(INTERNAL) Reset elliptic function and return M{e2 / f1 * cH * ...}.
        '''
        self._eF_reset_k2(x)
        return self._e2_f1 * y * self._eF.cH

    def _eF_reset_k2(self, x):
        '''(INTERNAL) Reset elliptic function and return C{k2}.
        '''
        ep2 = self.ep2
        k2 = x**2 * ep2
        self._eF.reset(k2=-k2,        alpha2=-ep2,
                       kp2=k2 + _1_0, alphap2=ep2 + _1_0)  # defaults
        _update_glXs(self)  # zap cached/memoized _GeodesicLineExact attrs
        return k2

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self._E

    @Property_RO
    def ep2(self):  # aka .e22
        '''Get the ellipsoid's I{2nd eccentricity squared} (C{float}), M{e2 / (1 - e2)}.
        '''
        return self.ellipsoid.e22  # == self.e2 / self.f1**2

    @Property_RO
    def _eTOL2(self):
        '''(INTERNAL) The si12 threshold for "really short".
        '''
        return _eTOL2(self.f)

    @Property_RO
    def f(self):
        '''Get the ellipsoid's I{flattening} (C{float}), M{(a - b) / a}, C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    flattening = f

    @Property_RO
    def f1(self):  # aka .b_a
        '''Get the ellipsoid's ratio I{polar} over I{equatorial} radius (C{float}), M{b / a} == M{1 - f}.
        '''
        return self.ellipsoid.b_a  # _1_0 - self.f

    @Property_RO
    def _f_180(self):
        '''(INTERNAL) Cached/memoized.
        '''
        return self.f * _180_0

    def _GDictDirect(self, lat1, lon1, azi1, arcmode, s12_a12, outmask=Caps.STANDARD):
        '''(INTERNAL) As C{_GenDirect}, but returning a L{GDict}.

           @return: A L{GDict} ...
        '''
        if not arcmode:  # supply DISTANCE_IN
            outmask |= Caps.DISTANCE_IN
        glX = self._LineTemp(lat1, lon1, azi1, outmask)
        return glX._GDictPosition(arcmode, s12_a12, outmask)

    def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=Caps.STANDARD):  # MCCABE 33, 41 vars
        '''(INTERNAL) As C{_GenInverse}, but returning a L{GDict}.

           @return: A L{GDict} ...
        '''
        if self.debug:  # PYCHOK no cover
            outmask |= Caps._DEBUG_INVERSE & self._debug
        outmask &= Caps._OUTMASK  # incl. _SALPs_CALPs and _DEBUG_
        # compute longitude difference carefully (with _diff182):
        # result is in [-180, +180] but -180 is only for west-going
        # geodesics, +180 is for east-going and meridional geodesics
        lon12, lon12s = _diff182(lon1, lon2)
        # see C{result} from geographiclib.geodesic.Inverse
        if (outmask & Caps.LONG_UNROLL):  # == (lon1 + lon12) + lon12s
            r = GDict(lon1=lon1, lon2=fsum_(lon1, lon12, lon12s))
        else:
            r = GDict(lon1=_norm180(lon1), lon2=_norm180(lon2))
        # make longitude difference positive and if very close to
        # being on the same half-meridian, then make it so.
        if lon12 < 0:
            lon_, lon12 = True, -_around(lon12)
            lon12s = _around((_180_0 - lon12) + lon12s)
        else:
            lon_, lon12 = False, _around(lon12)
            lon12s = _around((_180_0 - lon12) - lon12s)
        lam12 = radians(lon12)
        if lon12 > _90_0:
            slam12, clam12 = sincos2d(lon12s)
            clam12 = -clam12
        else:
            slam12, clam12 = sincos2d(lon12)

        # If really close to the equator, treat as on equator.
        lat1 = _around(_fix90(lat1))
        lat2 = _around(_fix90(lat2))
        r.set_(lat1=lat1, lat2=lat2)
        # Swap points so that point with higher (abs) latitude is
        # point 1.  If one latitude is a NAN, then it becomes lat1.
        swap_ = abs(lat1) < abs(lat2)
        if swap_:
            lat1, lat2 = lat2, lat1
            lon_ = not lon_
        if lat1 < 0:
            lat_ = False
        else:  # make lat1 <= 0
            lat_ = True
            lat1, lat2 = -lat1, -lat2
        # Now 0 <= lon12 <= 180, -90 <= lat1 <= 0 and lat1 <= lat2 <= -lat1
        # and lat_, lon_, swap_ register the transformation to bring the
        # coordinates to this canonical form, where False means no change
        # made.  We make these transformations so that there are few cases
        # to check, e.g., on verifying quadrants in atan2.  In addition,
        # this enforces some symmetries in the results returned.

        # Initialize for the meridian.  No longitude calculation is
        # done in this case to let the parameter default to 0.
        sbet1, cbet1 = self._sincos2f1(lat1)
        sbet2, cbet2 = self._sincos2f1(lat2)
        # If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure
        # of the |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1),
        # abs(sbet2) + sbet1 is a better measure.  This logic is used
        # in assigning calp2 in _Lambda10.  Sometimes these quantities
        # vanish and in that case we force bet2 = +/- bet1 exactly.  An
        # example where is is necessary is the inverse problem
        # 48.522876735459 0 -48.52287673545898293 179.599720456223079643
        # which failed with Visual Studio 10 (Release and Debug)
        if cbet1 < -sbet1:
            if cbet2 == cbet1:
                sbet2 = sbet1 if sbet2 < 0 else -sbet1
        elif abs(sbet2) == -sbet1:
            cbet2 = cbet1

        p = _PDict(sbet1=sbet1, cbet1=cbet1, dn1=self._dn(sbet1, cbet1),
                   sbet2=sbet2, cbet2=cbet2, dn2=self._dn(sbet2, cbet2))

        _meridian = True  # i.e. meridian = False
        if lat1 == -90 or slam12 == 0:
            # Endpoints are on a single full meridian,
            # so the geodesic might lie on a meridian.
            salp1, calp1 =  slam12, clam12  # head to target lon
            salp2, calp2 = _0_0,   _1_0  # then head north
            # tan(bet) = tan(sig) * cos(alp)
            p.setsigs(sbet1, calp1 * cbet1, sbet2, calp2 * cbet2)
            # sig12 = sig2 - sig1
            sig12 = atan2(*_sincos12(sbet1, p.csig1, sbet2, p.csig2, noneg=True))  # PYCHOK csig*
            s12x, m12x, _, \
            M12,  M21 = self._Lengths5(sig12, outmask | Caps.REDUCEDLENGTH, p)
            # Add the check for sig12 since zero length geodesics
            # might yield m12 < 0.  Test case was
            #    echo 20.001 0 20.001 0 | GeodSolve -i
            # In fact, we will have sig12 > PI/2 for meridional
            # geodesic which is not a shortest path.
            if m12x >= 0 or sig12 < _1_0:
                # Need at least 2 to handle 90 0 90 180
                # Prevent negative s12 or m12 from geographiclib 1.52
                if sig12 < _TINY3 or (sig12 < _TOL0 and (s12x < 0 or m12x < 0)):
                    sig12 = m12x = s12x = a12 = _0_0
                else:
                    m12x *= self.b
                    s12x *= self.b
                    a12 = degrees(sig12)
                _meridian = False
                C = 1
            # elif m12 < 0:  # prolate and too close to anti-podal
            #   _meridian = True

        if _meridian:
            if sbet1 == 0 and (  # sbet2 == 0 and
               self.f <= 0 or lon12s >= self._f_180):
                # Geodesic runs along equator
                calp1 = calp2 = _0_0
                salp1 = salp2 = _1_0
                sig12 = lam12 / self.f1  # == omg12
                somg12, comg12 = sincos2(sig12)
                m12x = self.b * somg12
                s12x = self.a * lam12
                if (outmask & Caps.GEODESICSCALE):
                    r.set_(M12=comg12, M21=comg12)
                a12 = lon12 / self.f1
                C = 2
            else:
                # Now point1 and point2 belong within a hemisphere bounded by a
                # meridian and geodesic is neither meridional or equatorial.
                p.set_(slam12=slam12, clam12=clam12)
                # Figure a starting point for Newton's method
                sig12, salp1, calp1, \
                       salp2, calp2, dnm = self._InverseStart6(lam12, p)
                if sig12 is None:  # use Newton's method
                    # pre-compute the constant _Lambda10 term, once
                    p.set_(bet12=None if cbet2 == cbet1 and abs(sbet2) == -sbet1 else (
                               ((cbet1 + cbet2) * (cbet2 - cbet1)) if cbet1 < -sbet1 else
                               ((sbet1 + sbet2) * (sbet1 - sbet2))))
                    sig12, salp1, calp1, salp2, calp2, domg12, \
                           s12x, m12x, M12, M21 = self._Newton10(salp1, calp1, outmask, p)
                    if (outmask & Caps.AREA):
                        # omg12 = lam12 - domg12
                        s, c = sincos2(domg12)
                        somg12, comg12 = _sincos12(s, c, slam12, clam12)

                    if (outmask & Caps.GEODESICSCALE):
                        if swap_:
                            M12, M21 = M21, M12
                        r.set_(M12=unsigned0(M12),
                               M21=unsigned0(M21))
                    C = 3  # Newton
                else:
                    C = 4  # Short lines, _Inverse6 set salp2, calp2, dnm
                    s, c = sincos2(sig12 / dnm)
                    m12x = s * dnm**2
                    s12x = sig12 * dnm
                    if (outmask & Caps.GEODESICSCALE):
                        r.set_(M12=c, M21=c)
                    if (outmask & Caps.AREA):
                        somg12, comg12 = sincos2(lam12 / (self.f1 * dnm))
                a12   = degrees(sig12)
                m12x *= self.b
                s12x *= self.b

        else:  # _meridian is False
            somg12 = comg12 = NAN

        r.set_(a12=a12)  # in [0, 180]
        if outmask == Caps._ANGLE_ONLY:
            return r  # for .Inverse1

        if (outmask & Caps.DISTANCE):
            r.set_(s12=unsigned0(s12x))

        if (outmask & Caps.REDUCEDLENGTH):
            r.set_(m12=unsigned0(m12x))

        if (outmask & Caps.AREA):
            S12 = self._InverseArea(_meridian, salp1, calp1,
                                               salp2, calp2,
                                               somg12, comg12, p)
            if _xor(swap_, lat_, lon_):
                S12 = -S12
            r.set_(S12=unsigned0(S12))

        if swap_:
            salp1, salp2 = salp2, salp1
            calp1, calp2 = calp2, calp1
        if _xor(swap_, lon_):
            salp1, salp2 = -salp1, -salp2
        if _xor(swap_, lat_):
            calp1, calp2 = -calp1, -calp2

        if (outmask & Caps.AZIMUTH):
            r.set_(azi1=atan2d(salp1, calp1),
                   azi2=atan2d(salp2, calp2, reverse=outmask & Caps.REVERSE2))
        if (outmask & Caps._SALPs_CALPs):
            r.set_(salp1=salp1, calp1=calp1,
                   salp2=salp2, calp2=calp2)

        if (outmask & Caps._DEBUG_INVERSE):  # PYCHOK no cover
            eF = self._eF
            p.set_(C=C, a=self.a, f=self.f, f1=self.f1,
                        e=self.ellipsoid.e, e2=self.e2, ep2=self.ep2,
                        c2x=self.c2x, c2=self.ellipsoid.c2,
                        eFcD=eF.cD, eFcE=eF.cE, eFcH=eF.cH,
                        eFk2=eF.k2, eFa2=eF.alpha2)
            p.update(r)  # r overrides
            r = p.toGDict()
        return r

    def _GenDirect(self, lat1, lon1, azi1, arcmode, s12_a12, outmask):
        '''(INTERNAL) The general I{Inverse} geodesic calculation.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2,
                                      s12, m12,  M12,  M21, S12)}.
        '''
        r = self._GDictDirect(lat1, lon1, azi1, arcmode, s12_a12, outmask)
        return r.toDirect9Tuple()

    def _GenDirectLine(self, lat1, lon1, azi1, arcmode, s12_a12, caps):
        '''(INTERNAL) Helper for C{ArcDirectLine} and C{DirectLine}.

           @return: A L{GeodesicLineExact} instance.
        '''
        azi1 = _norm180(azi1)
        # guard against underflow in salp0.  Also -0 is converted to +0.
        s, c = sincos2d(_around(azi1))

        if not arcmode:  # supply DISTANCE_IN
            caps |= Caps.DISTANCE_IN
        return _GeodesicLineExact(self, lat1, lon1, azi1, caps,
                                  self._debug, s, c)._GenSet(arcmode, s12_a12)

    def _GenInverse(self, lat1, lon1, lat2, lon2, outmask):
        '''(INTERNAL) The general I{Inverse} geodesic calculation.

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1, salp2, calp2,
                                             m12,   M12,   M21,   S12)}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, outmask | Caps._SALPs_CALPs)
        return r.toInverse10Tuple()

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.STANDARD):
        '''Perform the I{Inverse} geodesic calculation.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg lat2: Latitude of the second point (C{degrees}).
           @arg lon2: Longitude of the second point (C{degrees}).
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                    lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                    C{lon1}, C{azi1} and distance C{s12} always included.

           @note: The third point of the L{GeodesicLineExact} is set to correspond
                  to the second point of the I{Inverse} geodesic problem.

           @note: Both B{C{lat1}} and B{C{lat2}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.InverseLine
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.InverseLine<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return self._GDictInverse(lat1, lon1, lat2, lon2, outmask)

    def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
        '''Return the non-negative, I{angular} distance in C{degrees}.
        '''
        # see .FrechetKarney.distance, .HausdorffKarney._distance
        # and .HeightIDWkarney._distances
        _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
        return abs(self._GDictInverse(lat1, lon1, lat2, lon2, Caps._ANGLE_ONLY).a12)

    def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
        '''Return the distance in C{meter} and the forward and
           reverse azimuths (initial and final bearing) in C{degrees}.

           @return: L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, Caps._AZIMUTH_DISTANCE)
        return Distance3Tuple(r.s12, wrap360(r.azi1), wrap360(r.azi2))

    def InverseLine(self, lat1, lon1, lat2, lon2, caps=Caps.STANDARD):
        '''Define a L{GeodesicLineExact} in terms of the I{Inverse} geodesic problem.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg lat2: Latitude of the second point (C{degrees}).
           @arg lon2: Longitude of the second point (C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineExact} instance
                        should possess, i.e., which quantities can be
                        returned by calls to L{GeodesicLineExact.Position}
                        and L{GeodesicLineExact.ArcPosition}.

           @return: A L{GeodesicLineExact} instance.

           @note: The third point of the L{GeodesicLineExact} is set to correspond
                  to the second point of the I{Inverse} geodesic problem.

           @note: Both B{C{lat1}} and B{C{lat2}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.InverseLine
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.InverseLine<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, Caps._SALPs_CALPs)  # No need for AZIMUTH
        azi1 = atan2d(r.salp1, r.calp1)
        if (caps & (Caps._OUTMASK & Caps.DISTANCE_IN)):
            caps |= Caps.DISTANCE  # ensure that a12 is a distance
        return _GeodesicLineExact(self, lat1, lon1, azi1, caps,
                                  self._debug, r.salp1, r.calp1)._GenSet(True, r.a12)

    def _InverseArea(self, _meridian, salp1, calp1,  # PYCHOK 9 args
                                      salp2, calp2,
                                      somg12, comg12, p):
        '''(INTERNAL) Split off from C{_GDictInverse} to reduce complexity/length.

           @return: Area C{S12}.
        '''
        # from _Lambda10: sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * p.cbet1
        calp0 = hypot(calp1, salp1 * p.sbet1)  # calp0 > 0
        if calp0 and salp0:
            # from _Lambda10: tan(bet) = tan(sig) * cos(alp)
            k2  =  calp0**2 * self.ep2
            C4a =  self._C4f_k2(k2)
            B41 = _coSeries(C4a, *norm2(p.sbet1, calp1 * p.cbet1))
            B42 = _coSeries(C4a, *norm2(p.sbet2, calp2 * p.cbet2))
            # multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
            A4  = self._e2a2 * calp0 * salp0
            S12 = A4 * (B42 - B41)
            p.set_(A4=A4, B41=B41, B42=B42, k2=k2)
        else:  # avoid problems with indeterminate sig1, sig2 on equator
            S12 = _0_0

        if (_meridian and  # omg12 < 3/4 * PI
             comg12 > _SQRT2_2 and  # lon diff not too big
             (p.sbet2 - p.sbet1) < _1_75):  # lat diff not too big
            # use tan(Gamma/2) = tan(omg12/2) *
            #                   (tan(bet1/2) + tan(bet2/2)) /
            #                   (tan(bet1/2) * tan(bet2/2) + 1))
            # with tan(x/2) = sin(x) / (1 + cos(x))
            cbet1  = _1_0 + p.cbet1
            cbet2  = _1_0 + p.cbet2
            salp12 = (p.sbet1 *   cbet2 + cbet1 * p.sbet2) *  somg12
            calp12 = (p.sbet1 * p.sbet2 + cbet1 *   cbet2) * (comg12 + _1_0)
            alp12  = _2_0 * atan2(salp12, calp12)
        else:
            # alp12 = alp2 - alp1, used in atan2, no need to normalize
            salp12, calp12 = _sincos12(salp1, calp1, salp2, calp2)
            # The right thing appears to happen if alp1 = +/-180 and
            # alp2 = 0, viz salp12 = -0 and alp12 = -180.  However,
            # this depends on the sign being attached to 0 correctly.
            # Following ensures the correct behavior.
            if salp12 == 0 and calp12 < 0:
                alp12 = copysign0(PI, calp1)
            else:
                alp12 = atan2(salp12, calp12)

        p.set_(alp12=alp12)
        return S12 + self.c2x * alp12

    def _InverseStart6(self, lam12, p):
        '''(INTERNAL) Return a starting point for Newton's method in
           C{salp1} and C{calp1} indicated by C{sig12=None}.  If
           Newton's method doesn't need to be used, return also
           C{salp2}, C{calp2}, C{dnm} and C{sig12} non-C{None}.

           @return: 6-Tuple C{(sig12, salp1, calp1, salp2, calp2, dnm)}
                    and updated C{p.setsigs} if C{sig12==None}.
        '''
        sig12 = None  # use Newton
        salp1 = calp1 = salp2 = calp2 = dnm = NAN

        # bet12 = bet2 - bet1 in [0, PI)
        sbet12, cbet12 = _sincos12(p.sbet1, p.cbet1, p.sbet2, p.cbet2)
        shortline = cbet12 >= 0 and sbet12 < _0_5 and (p.cbet2 * lam12) < _0_5
        if shortline:
            # sin((bet1 + bet2)/2)^2 = (sbet1 + sbet2)^2 / (
            #      (cbet1 + cbet2)^2 + (sbet1 + sbet2)^2)
            t = (p.sbet1 + p.sbet2)**2
            s = t / ((p.cbet1 + p.cbet2)**2 + t)
            dnm = sqrt(_1_0 + self.ep2 * s)
            somg12, comg12 = sincos2(lam12 / (self.f1 * dnm))
        else:
            somg12, comg12 = p.slam12, p.clam12

        # bet12a = bet2 + bet1 in (-PI, 0], note -sbet1
        sbet12a, cbet12a = _sincos12(-p.sbet1, p.cbet1, p.sbet2, p.cbet2)

        s = somg12**2 / (abs(comg12) + _1_0)
        t = p.cbet2 * p.sbet1 * s
        salp1 =  p.cbet2 * somg12
        calp1 = (sbet12a - t) if comg12 < 0 else (sbet12 + t)

        ssig12 = hypot(salp1, calp1)
        csig12 = p.sbet1 * p.sbet2 + p.cbet1 * p.cbet2 * comg12

        if shortline and ssig12 < self._eTOL2:  # really short lines
            if comg12 < 0:
                s = _1_0 - comg12
            salp2, calp2 = norm2(somg12 * p.cbet1,
                                 sbet12 - p.cbet1 * p.sbet2 * s)
            sig12 = atan2(ssig12, csig12)  # return value

        elif (self._n_0_1 or  # Skip astroid calc if too eccentric
              csig12 >= 0 or ssig12 >= (self._n_6PI * p.cbet1**2)):
            pass  # nothing to do, 0th order spherical approximation OK

        else:
            # Scale lam12 and bet2 to x, y coordinate system where antipodal
            # point is at origin and singular point is at y = 0, x = -1
            lam12x = atan2(-p.slam12, -p.clam12)  # lam12 - PI
            f = self.f
            if f < 0:  # x = dlat, y = dlon
                # ssig1=sbet1, csig1=-cbet1, ssig2=sbet2, csig2=cbet2
                p.setsigs(p.sbet1, -p.cbet1, p.sbet2, p.cbet2)
                # if lon12 = 180, this repeats a calculation made in Inverse
                _, m12b, m0, _, _ = self._Lengths5(atan2(sbet12a, cbet12a) + PI,
                                                   Caps.REDUCEDLENGTH, p)
                t = p.cbet1 * PI
                x = m12b / (t * p.cbet2 * m0) - _1_0
                sca = (sbet12a / (x * p.cbet1)) if x < -_0_01 else (-f * t)
                y = lam12x / sca
            else:  # _f >= 0, in fact f == 0 does not get here
                sca = self._eF_reset_e2_f1_cH(p.sbet1, p.cbet1 * _2_0)
                x = lam12x / sca  # dlon
                y = sbet12a / (sca * p.cbet1)  # dlat

            if y > _TOL1 and x > -_THR1:  # strip near cut
                if f < 0:
                    calp1 = max((_0_0 if x > _TOL1 else _N_1_0), x)
                    salp1 = sqrt(_1_0 - calp1**2)
                else:
                    salp1 =  min(_1_0, -x)
                    calp1 = -sqrt(_1_0 - salp1**2)
            else:
                # Estimate alp1, by solving the astroid problem.
                #
                # Could estimate alpha1 = theta + PI/2, directly, i.e.,
                #   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
                #   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
                #
                # However, it's better to estimate omg12 from astroid and use
                # spherical formula to compute alp1.  This reduces the mean
                # number of Newton iterations for astroid cases from 2.24
                # (min 0, max 6) to 2.12 (min 0 max 5).  The changes in the
                # number of iterations are as follows:
                #
                # change percent
                #    1       5
                #    0      78
                #   -1      16
                #   -2       0.6
                #   -3       0.04
                #   -4       0.002
                #
                # The histogram of iterations is (m = number of iterations
                # estimating alp1 directly, n = number of iterations
                # estimating via omg12, total number of trials = 148605):
                #
                #  iter    m      n
                #    0   148    186
                #    1 13046  13845
                #    2 93315 102225
                #    3 36189  32341
                #    4  5396      7
                #    5   455      1
                #    6    56      0
                #
                # omg12 is near PI, estimate work with omg12a = PI - omg12
                k = _Astroid(x, y)
                k1 = k + _1_0
                sca *= (y * k1 / k) if f < 0 else (x * k / k1)
                s, c = sincos2(-sca)
                # update spherical estimate of alp1 using omg12 instead of lam12
                salp1 = p.cbet2 * s
                calp1 = sbet12a - p.cbet2 * p.sbet1 * s**2 / (_1_0 + c)

        # sanity check on starting guess.  Backwards check allows NaN through.
        salp1, calp1 = norm2(salp1, calp1) if salp1 > 0 else (_1_0, _0_0)

        return sig12, salp1, calp1, salp2, calp2, dnm

    def _Lambda5(self, salp1, calp1, diffp, p):
        '''(INTERNAL) Helper.

           @return: 5-Tuple C{(lam12, sig12, salp2, calp2, dlam12)}
                    and updated C{p.setsigs} and C{p.domg12}.
        '''
        eF = self._eF
        f1 = self.f1

        if p.sbet1 == 0 and calp1 == 0:
            # Break degeneracy of equatorial line
            calp1 = -_TINY

        # sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * p.cbet1
        calp0 = hypot(calp1, salp1 * p.sbet1)  # calp0 > 0
        # tan(bet1) = tan(sig1) * cos(alp1)
        # tan(omg1) = sin(alp0) * tan(sig1)
        #           = sin(bet1) * tan(alp1)
        somg1 = salp0 * p.sbet1
        comg1 = calp1 * p.cbet1
        ssig1, csig1 = norm2(p.sbet1, comg1)
        # Without normalization we have schi1 = somg1
        cchi1 = f1 * p.dn1 * comg1

        # Enforce symmetries in the case abs(bet2) = -bet1.
        # Need to be careful about this case, since this can
        # yield singularities in the Newton iteration.
        # sin(alp2) * cos(bet2) = sin(alp0)
        salp2 = (salp0 / p.cbet2) if p.cbet2 != p.cbet1 else salp1
        # calp2 = sqrt(1 - sq(salp2))
        #       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
        # and subst for calp0 and rearrange to give (choose
        # positive sqrt to give alp2 in [0, PI/2]).
        calp2 = abs(calp1) if p.bet12 is None else (
                sqrt(p.bet12 + (calp1 * p.cbet1)**2) / p.cbet2)
        # tan(bet2) = tan(sig2) * cos(alp2)
        # tan(omg2) = sin(alp0) * tan(sig2).
        somg2 = salp0 * p.sbet2
        comg2 = calp2 * p.cbet2
        ssig2, csig2 = norm2(p.sbet2, comg2)
        # without normalization we have schi2 = somg2
        cchi2 = f1 * p.dn2 * comg2

        # omg12 = omg2 - omg1, limit to [0, PI]
        somg12, comg12 = _sincos12(somg1, comg1, somg2, comg2, noneg=True)
        # chi12 = chi2 - chi1, limit to [0, PI]
        schi12, cchi12 = _sincos12(somg1, cchi1, somg2, cchi2, noneg=True)

        p.setsigs(ssig1, csig1, ssig2, csig2)
        # sig12 = sig2 - sig1, limit to [0, PI]
        sig12  = atan2(*_sincos12(ssig1, csig1, ssig2, csig2, noneg=True))

        eta12  = self._eF_reset_e2_f1_cH(calp0, salp0) / PI_2  # then ...
        eta12 *= fsum_(eF.deltaH(*p.sncndn2),
                      -eF.deltaH(*p.sncndn1), sig12)
        # eta = chi12 - lam12
        lam12  = atan2(*_sincos12(p.slam12, p.clam12, schi12, cchi12)) - eta12
        # domg12 = chi12 - omg12 - deta12
        domg12 = atan2(*_sincos12(  somg12,   comg12, schi12, cchi12)) - eta12

        if not diffp:
            dlam12  = NAN  # dv > 0 in ._Newton9
        elif calp2:
            _, dlam12, _, _, _ = self._Lengths5(sig12, Caps.REDUCEDLENGTH, p)
            dlam12 *=  f1 / (calp2 * p.cbet2)
        else:
            dlam12  = -f1 * _2_0 * p.dn1 / p.sbet1

        # p.set_(deta12=-eta12, lam12=lam12)
        return lam12, sig12, salp2, calp2, domg12, dlam12

    def _Lengths5(self, sig12, outmask, p):
        '''(INTERNAL) Return M{m12b = (reduced length) / self.b} and
           calculate M{s12b = distance / self.b} and M{m0}, the
           coefficient of secular term in expression for reduced length.

           @return: 5-Tuple C{(s12b, m12b, m0, M12, M21)}.
        '''
        s12b = m12b = m0 = M12 = M21 = NAN

        eF = self._eF

        # outmask &= Caps._OUTMASK
        if (outmask & Caps.DISTANCE):
            # Missing a factor of self.b
            s12b = eF.cE / PI_2 * fsum_(eF.deltaE(*p.sncndn2),
                                       -eF.deltaE(*p.sncndn1), sig12)

        if (outmask & Caps._REDUCEDLENGTH_GEODESICSCALE):
            m0x = -eF.k2 * eF.cD / PI_2
            J12 = -m0x * fsum_(eF.deltaD(*p.sncndn2),
                              -eF.deltaD(*p.sncndn1), sig12)
            c12 = p.csig1 * p.csig2
            if (outmask & Caps.REDUCEDLENGTH):
                m0 = m0x
                # Missing a factor of self.b.  Add parens around
                # (csig1 * ssig2) and (ssig1 * csig2) to ensure
                # accurate cancellation for coincident points.
                m12b = fsum_(p.dn2 * (p.csig1 * p.ssig2),
                            -p.dn1 * (p.ssig1 * p.csig2), J12 * c12)

            if (outmask & Caps.GEODESICSCALE):
                M12 = M21 = p.ssig1 * p.ssig2 + c12
                t = (p.cbet1 - p.cbet2) * self.ep2 * \
                    (p.cbet1 + p.cbet2) / (p.dn1 + p.dn2)
                M12 += (p.ssig2 * t + p.csig2 * J12) * p.ssig1 / p.dn1
                M21 -= (p.ssig1 * t + p.csig1 * J12) * p.ssig2 / p.dn2

        return s12b, m12b, m0, M12, M21

    def Line(self, lat1, lon1, azi1, caps=Caps.ALL):
        '''Set up a L{GeodesicLineExact} to compute several points
           on a single geodesic.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineExact} instance
                        should possess, i.e., which quantities can be
                        returnedby calls to L{GeodesicLineExact.Position}
                        and L{GeodesicLineExact.ArcPosition}.

           @return: A L{GeodesicLineExact} instance.

           @note: If the point is at a pole, the azimuth is defined by keeping
                  B{C{lon1}} fixed, writing C{B{lat1} = ±(90 − ε)}, and taking
                  the limit C{ε → 0+}.

           @see: C++ U{GeodesicExact.Line
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Line<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return _GeodesicLineExact(self, lat1, lon1, azi1, caps, self._debug)

    def _LineTemp(self, lat1, lon1, azi1, caps):  # for ._GDictDirect and in .azimuthal._GnomonicBase.reverse
        '''(INTERNAL) Get a temporary L{GeodesicLineExact} instance.
        '''
        return _GeodesicLineExact(self, lat1, lon1, azi1, caps, self._debug, append=False)

    @Property_RO
    def n(self):
        '''Get the ellipsoid's I{3rd flattening} (C{float}), M{f / (2 - f) == (a - b) / (a + b)}.
        '''
        return self.ellipsoid.n

    @Property_RO
    def _n_0_1(self):
        '''(INTERNAL) Cached once.
        '''
        return abs(self.n) > _0_1

    @Property_RO
    def _n_6PI(self):
        '''(INTERNAL) Cached once.
        '''
        return abs(self.n) * _6_0 * PI

    def _Newton10(self, salp1, calp1, outmask, p):
        '''(INTERNAL) Split off from C{_GDictInverse} to reduce complexity/length.

           @return: 10-Tuple C{(sig12, salp1, calp1, salp2, calp2,
                    domg12, s12x, m12x, M12, M21)}.
        '''
        # This is a straightforward solution of f(alp1) = lambda12(alp1) -
        # lam12 = 0 with one wrinkle.  f(alp) has exactly one root in the
        # interval (0, PI) and its derivative is positive at the root.
        # Thus f(alp) is positive for alp > alp1 and negative for alp < alp1.
        # During the course of the iteration, a range (alp1a, alp1b) is
        # maintained which brackets the root and with each evaluation of
        # f(alp) the range is shrunk, if possible.  Newton's method is
        # restarted whenever the derivative of f is negative (because the
        # new value of alp1 is then further from the solution) or if the
        # new estimate of alp1 lies outside (0,PI); in this case, the new
        # starting guess is taken to be (alp1a + alp1b) / 2.
        salp1a = salp1b = _TINY
        calp1a,  calp1b = _1_0, _N_1_0
        tripb,   TOLv   =  False, _TOL0
        for it in range(_MAXIT2):
            # 1/4 meridian = 10e6 meter and random input,
            # estimated max error in nm (nano meter, by
            # checking Inverse problem by Direct).
            #
            #             max   iterations
            # log2(b/a)  error  mean   sd
            #    -7       387   5.33  3.68
            #    -6       345   5.19  3.43
            #    -5       269   5.00  3.05
            #    -4       210   4.76  2.44
            #    -3       115   4.55  1.87
            #    -2        69   4.35  1.38
            #    -1        36   4.05  1.03
            #     0        15   0.01  0.13
            #     1        25   5.10  1.53
            #     2        96   5.61  2.09
            #     3       318   6.02  2.74
            #     4       985   6.24  3.22
            #     5      2352   6.32  3.44
            #     6      6008   6.30  3.45
            #     7     19024   6.19  3.30
            v, sig12, salp2, calp2, \
               domg12, dv = self._Lambda5(salp1, calp1, it < _MAXIT1, p)

            # 2 * _TOL0 is approximately 1 ulp [0, PI]
            # reversed test to allow escape with NaNs
            if tripb or abs(v) < TOLv:
                break
            # update bracketing values
            if v > 0 and (it > _MAXIT1 or (calp1 / salp1) > (calp1b / salp1b)):
                salp1b, calp1b = salp1, calp1
            elif v < 0 and (it > _MAXIT1 or (calp1 / salp1) < (calp1a / salp1a)):
                salp1a, calp1a = salp1, calp1

            if it < _MAXIT1 and dv > 0:
                dalp1 = -v / dv
                if abs(dalp1) < PI:
                    s, c = sincos2(dalp1)
                    # nalp1 = alp1 + dalp1
                    s, c = _sincos12(-s, c, salp1, calp1)
                    if s > 0:
                        salp1, calp1 = norm2(s, c)
                        # in some regimes we don't get quadratic convergence
                        # because slope -> 0.  So use convergence conditions
                        # based on epsilon instead of sqrt(epsilon)
                        TOLv = _TOL0 if abs(v) > _TOL016 else _TOL08
                        continue

            # Either dv was not positive or updated value was outside
            # legal range.  Use the midpoint of the bracket as the next
            # estimate.  This mechanism is not needed for the WGS84
            # ellipsoid, but it does catch problems with more eccentric
            # ellipsoids.  Its efficacy is such for the WGS84 test set
            # with the starting guess set to alp1 = 90deg: the WGS84
            # test set: mean = 5.21, sd = 3.93, max = 24 WGS84 and
            # random input: mean = 4.74, sd = 0.99
            salp1, calp1 = norm2((salp1a + salp1b) * _0_5,
                                 (calp1a + calp1b) * _0_5)
            tripb = fsum_(calp1a, -calp1, abs(salp1a - salp1)) < _TOLb or \
                    fsum_(calp1b, -calp1, abs(salp1b - salp1)) < _TOLb
            TOLv = _TOL0

        else:
            raise GeodesicError(_no_(_convergence_), _MAXIT2, txt=repr(self))  # self.toRepr()

        s12x, m12x, _, M12, M21 = self._Lengths5(sig12, outmask, p)

        p.set_(iter=it, trip=tripb)
        return (sig12, salp1, calp1, salp2, calp2, domg12,
                       s12x,  m12x,  M12,   M21)

    def _sincos2f1(self, lat):
        '''(INTERNAL) Helper.
        '''
        sbet, cbet = sincos2d(lat)
        # ensure cbet1 = +epsilon at poles; doing the fix on beta means
        # that sig12 will be <= 2*tiny for two points at the same pole
        sbet, cbet = norm2(sbet * self.f1, cbet)
        return sbet, max(_TINY, cbet)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{GeodesicExact} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).

           @return: Tuple items (C{str}).
        '''
        d = dict(ellipsoid=self.ellipsoid, C4Order=self.C4Order)
        return sep.join(pairs(d, prec=prec))


class GeodesicLineExact(_GeodesicLineExact):
    '''A pure Python version of I{Karney}'s C++ class U{GeodesicLineExact
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicLineExact.html>},
       modeled after I{Karney}'s Python class U{GeodesicLine<https://GeographicLib.SourceForge.io/
       html/python/code.html#module-geographiclib.geodesicline>}.
    '''

    def __init__(self, geodesic, lat1, lon1, azi1, caps=Caps.STANDARD, name=NN):
        '''New L{GeodesicLineExact} instance, allowing points to be found along
           a geodesic starting at C{(B{lat1}, B{lon1})} with azimuth B{C{azi1}}.

           @arg geodesic: The geodesic to use (L{GeodesicExact}).
           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first points (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineExact} instance
                        should possess, i.e., which quantities can be
                        returned by calls to L{GeodesicLineExact.Position}
                        and L{GeodesicLineExact.ArcPosition}.
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{geodesic}}.
        '''
        _xinstanceof(GeodesicExact, geodesic=geodesic)

        _GeodesicLineExact.__init__(self, geodesic, lat1, lon1, azi1, caps, 0, name=name)


def _Astroid(x, y):
    '''(INTERNAL) Solve M{k^4 + 2 * k^3 - (x^2 + y^2 - 1 ) * k^2 -
       (2 * k + 1) * y^2 = 0} for positive root k.
    '''
    p = x**2
    q = y**2
    r = q + p - _1_0
    if q or r > 0:
        r = r / _6_0  # /= chokes PyChecker
        # avoid possible division by zero when r = 0
        # by multiplying s and t by r^3 and r, resp.
        S = p * q / _4_0  # S = r^3 * s
        r3 = r**3
        T3 = r3 + S
        # discriminant of the quadratic equation for T3 is
        # zero on the evolute curve p^(1/3)+q^(1/3) = 1
        d = S * (S + _2_0 * r3)
        if d < 0:
            # T is complex, but u is defined for a real result
            a = atan2(sqrt(-d), -T3)
            # There are three possible cube roots, choose the one which
            # avoids cancellation.  Note d < 0 implies that r < 0.
            u = r * (_1_0 + _2_0 * cos(a / _3_0))
        else:
            # pick the sign on the sqrt to maximize abs(T3) to
            # minimize loss of precision due to cancellation.
            if d:
                T3 += copysign0(sqrt(d), T3)  # T3 = (r * t)^3
            # cbrt always returns the real root, cbrt(-8) = -2
            u = cbrt(T3)  # T = r * t
            if u:  # T can be zero; but then r2 / T -> 0
                u += r**2 / u
            u += r

        v = sqrt(u**2 + q)
        # avoid loss of accuracy when u < 0
        uv = (q / (v - u)) if u < 0 else (u + v)
        w = (uv - q) / (_2_0 * v)  # positive?
        # rearrange expression for k to avoid loss of accuracy due to
        # subtraction, division by 0 impossible because uv > 0, w >= 0
        k = uv / (sqrt(w**2 + uv) + w)  # guaranteed positive

    else:  # q == 0 && r <= 0
        # y = 0 with |x| <= 1.  Handle this case directly, for
        # y small, positive root is k = abs(y) / sqrt(1 - x^2)
        k = _0_0
    return k


__all__ += _ALL_DOCS(GeodesicExact, GeodesicLineExact)

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
