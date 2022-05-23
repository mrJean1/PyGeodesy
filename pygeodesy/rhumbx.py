
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ classes U{Rhumb
<https://GeographicLib.SourceForge.io/C++/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
<https://GeographicLib.SourceForge.io/C++/classGeographicLib_1_1RhumbLine.html>}.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, especially the U{Class List<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>},
the background information on U{Rhumb lines<https://GeographicLib.SourceForge.io/C++/doc/rhumb.html>},
the utily U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} and U{Online
rhumb line calculations<https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>}.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2014-2022)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isnan, istuplist, _xinstanceof
from pygeodesy.errors import _or, _ValueError, _xkwds
from pygeodesy.fmath import hypot, hypot1
from pygeodesy.fsums import fsum1_, pairs
from pygeodesy.interns import NAN, NN, PI_2, _COMMASPACE_, _lat1_, _lat2_, \
                             _lon1_, _lon2_, _not_, _s12_, _S12_, _UNDER_, \
                             _0_0, _0_5, _1_0, _2_0, _4_0, _90_0, _180_0, \
                             _720_0  # PYCHOK used!
from pygeodesy.karney import Caps, _CapsBase, _atan2d, _diff182, _ellipsoid, \
                             Direct9Tuple, _EWGS84, _fix90, GDict, _GTuple, \
                             Inverse10Tuple, _norm180, _polynomial, _tand, \
                             LatLon2Tuple
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
# from pygeodesy.namedTuples import LatLon2Tuple  # from .karney
from pygeodesy.props import Property, Property_RO, property_RO, _update_all
# from pygeodesy.streprs import pairs  # from .fsums
from pygeodesy.units import Bearing as _Azi, Int, Lat, Lon, \
                            Meter as _M, Meter2 as _M2
from pygeodesy.utily import sincos2_, sincos2d

from math import asinh, atan, cos, cosh, degrees, fabs, radians, sin, sinh, sqrt

__all__ = _ALL_LAZY.rhumbx + _ALL_DOCS(Caps)
__version__ = '22.05.22'

_azi12_ = 'azi12'

_rls = []  # instances of C{RbumbLine} to be updated


def _Lat(*lat, **Error_name):
    '''(INTERNAL) Latitude B{C{lat}}.
    '''
    kwds = _xkwds(Error_name, clip=0, Error=RhumbError)
    return Lat(*lat, **kwds)


def _Lon(*lon, **Error_name):
    '''(INTERNAL) Longitude B{C{lon}}.
    '''
    kwds = _xkwds(Error_name, clip=0, Error=RhumbError)
    return Lon(*lon, **kwds)


def _update_all_rls(r):
    '''(INTERNAL) Zap cached/memoized C{Property[_RO]}s
       of any L{RhumbLine} instances tied to the given
       L{Rhumb} instance B{C{r}}.
    '''
    _xinstanceof(r, Rhumb)
    _update_all(r)
    for rl in _rls:  # PYCHOK use weakref?
        if rl._rhumb is r:
            _update_all(rl)


class Rhumb(_CapsBase):
    '''Class to solve of the I{direct} and I{inverse rhumb} problems, accurately.

       @see: The U{Detailed Description<https://GeographicLib.SourceForge.io/C++/doc/
       classGeographicLib_1_1Rhumb.html>} of I{Karney}'s C++ C{Rhumb Class}.
    '''
    _E     = _EWGS84
    _exact =  True
    _mRA   =  6
    _mTM   =  6

    def __init__(self, a_earth=_EWGS84, f=None, exact=True, name=NN, **orders):
        '''New L{Rhumb}.

           @kwarg a_earth: This rhumb's earth (L{Ellipsoid}, L{Ellipsoid2},
                           L{a_f2Tuple}, L{Datum}, 2-tuple (C{a, f})) or the
                           equatorial radius (C{scalar}).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     a C{scalar}, ignored otherwise.
           @kwarg exact: If C{True}, use an addition theorem for elliptic integrals
                         to compute I{Divided differences}, otherwise use the Krüger
                         series expansion (C{bool}), see also property C{exact}.
           @kwarg name: Optional name (C{str}).
           @kwarg orders: Optional keyword arguments B{C{RAorder}} and B{C{TMorder}}
                          to set the respective C{order}, see properties C{RAorder}
                          and C{TMorder} and method C{orders}.

           @raise RhumbError: Invalid B{C{a_earth}}, B{C{f}} or B{C{orders}}.
        '''
        if f is not None:
            self.ellipsoid = a_earth, f
        elif a_earth not in (_EWGS84, None):
            self.ellipsoid = a_earth
        if not exact:
            self._exact = False
        if name:
            self.name = name
        if orders:
            self.orders(**orders)

    def _DConformal2Rectifying(self, chix, chiy):
        return _1_0 + _sincosSeries(True, chix, chiy,
                *self._DRectifying2ConformalAc2)

    @Property_RO
    def _DRectifying2ConformalAc2(self):
        return self._Xs2(_AlpCoeffs, self.TMorder)

    def Direct(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Solve the I{direct rhumb} problem, optionally with the area.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees}).
           @arg s12: Distance along the rhumb line from points 1 to
                     the second point (C{meter}), can be negative.

           @return: L{GDict} with 2 up to 9 items C{azi12, s12, S12, lat2,
                    lon2, lat1, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the other point in C{degrees} and the area
                    C{S12} under the rhumb line in C{meter} I{squared}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses
                  a pole, the longitude of the second point is indeterminate
                  and C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{epsilon}**-2 (where C{epsilon} is 2.0**-52.
                  This position is extremely close to the actual pole and
                  allows the calculation to be carried out in finite terms.
        '''
        rl = _RhumbLine(self, lat1, lon1, azi12, caps=Caps.LINE_OFF,
                              name=self.name)
        return rl.Position(s12, outmask | self._debug)  # lat2, lon2, S12

    def Direct7(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):
        '''Like method L{Rhumb.Direct} but returning a L{Rhumb7Tuple} with area C{S12}.
        '''
        return self.Direct(lat1, lon1, azi12, s12, outmask=outmask).toRhumb7Tuple()

    def DirectLine(self, lat1, lon1, azi12, name=NN):  # caps=Caps.STANDARD
        '''Define a L{RhumbLine} in terms of the I{direct} rhumb problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees}).

           @return: A L{RhumbLine} instance and invoke its method
                    L{RhumbLine.Position} to compute each point.

           @note: Updates to this rhumb are reflected in the returned
                  rhumb line.
        '''
        return RhumbLine(self, lat1=lat1, lon1=lon1, azi12=azi12,
                               caps=self._debug, name=name or self.name)

    def _DIsometric(self, latx, laty):
        phix = radians(latx)
        phiy = radians(laty)
        return (_Dtand(latx, laty) * _Dasinh(_tand(latx), _tand(laty))  -
                 _Dsin(phix, phiy) * _DeatanhE(sin(phix), sin(phiy), self._E))

    def _DIsometric2Rectifying(self, psix, psiy):  # degrees
        if self.exact:
            E = self._E
            latx = E.auxIsometric(psix, inverse=True)
            laty = E.auxIsometric(psiy, inverse=True)
            r = self._DRectifying(latx, laty) / self._DIsometric(latx, laty)
        else:
            x, y = radians(psix), radians(psiy)
            r = self._DConformal2Rectifying(_gd(x), _gd(y)) * _Dgd(x, y)
        return r

    def _DRectifying(self, latx, laty):
        E, eF = self._E, self._eF
        tbetx = E.f1 * _tand(latx)
        tbety = E.f1 * _tand(laty)
        return (E.f1 * _DfEt(tbetx, tbety, eF) * E.b
                     * _Dtand(latx, laty) * PI_2
                     * _Datan(tbetx, tbety)) / E.L

    def _DRectifying2Conformal(self, mux, muy):
        return _1_0 - _sincosSeries(True, mux, muy,
                *self._DRectifying2ConformalBc2)

    @Property_RO
    def _DRectifying2ConformalBc2(self):
        return self._Xs2(_BetCoeffs, self.TMorder)

    def _DRectifying2Isometric(self, mux, muy):  # radians
        E = self._E
        latx = E.auxRectifying(degrees(mux), inverse=True)
        laty = E.auxRectifying(degrees(muy), inverse=True)
        if self.exact:
            r = self._DIsometric(latx, laty) / self._DRectifying(latx, laty)
        else:
            r = self._DRectifying2Conformal(mux, muy) * \
                     _Dgdinv(E.es_taupf(_tand(latx)),
                             E.es_taupf(_tand(laty)))
        return r

    @Property_RO
    def _eF(self):
        '''(INTERNAL) Get the ellipsoid's elliptic function.
        '''
        # .k2 = 0.006739496742276434
        return self._E._elliptic_e12  # _MODS.elliptic.Elliptic(-self._E._e12)

    @Property
    def ellipsoid(self):
        '''Get this rhumb's ellipsoid (L{Ellipsoid}).
        '''
        return self._E

    @ellipsoid.setter  # PYCHOK setter!
    def ellipsoid(self, a_earth_f):
        '''Set this rhumb's ellipsoid (L{Ellipsoid}, L{Ellipsoid2}, L{Datum},
           L{a_f2Tuple} or 2-tuple C{(a, f)}).
        '''
        try:
            E = _ellipsoid(*a_earth_f[:2]) if istuplist(a_earth_f, 2) else \
                _ellipsoid( a_earth_f, None)
        except Exception as x:
            raise RhumbError(str(x))
        if E != self._E:
            _update_all_rls(self)
            self._E = E

    @Property
    def exact(self):
        '''Get the I{exact} option (C{bool}).
        '''
        return self._exact

    @exact.setter  # PYCHOK setter!
    def exact(self, exact):
        '''Set the I{exact} option (C{bool}).

           @arg exact: If C{True}, use I{exact} rhumb calculations, if C{False}
                       results are less precise for more oblate or more prolate
                       ellipsoids with M{abs(flattening) > 0.01} (C{bool}).

           @see: Option U{B{-s}<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}
                 and U{ACCURACY<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html#ACCURACY>}.
        '''
        x = bool(exact)
        if x != self.exact:
            _update_all_rls(self)
            self._exact = x

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH_DISTANCE):
        '''Solve the I{inverse rhumb} problem, optionally with the area.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).

           @return: L{GDict} with 4 to 7 itens C{lat1, lon2, lat2, lon2,
                    azi12, s12, S12}, the azimuth C{azi12} of the rhumb
                    line in compass C{degrees} between C{-180} and C{+180},
                    the rhumb C{s12} distance between both points in C{meter}
                    and the area C{S12} under the rhumb line in C{meter}
                    I{squared}.

           @note: The shortest rhumb line is found.  If the end points are
                  on opposite meridians, there are two shortest rhumb lines
                  and the East-going one is chosen.

                  If either point is a pole, the cosine of its latitude is
                  taken to be C{epsilon}**-2 (where C{epsilon} is 2.0**-52).
                  This position is extremely close to the actual pole and
                  allows the calculation to be carried out in finite terms.
        '''
        r = GDict(S12=NAN, name=self.name)
        if (outmask & Caps.AZIMUTH_DISTANCE_AREA):
            r.set_(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2)
            E = self._E
            psi1 = E.auxIsometric(lat1)
            psi2 = E.auxIsometric(lat2)
            psi12 = psi2 - psi1
            lon12, _ = _diff182(lon1, lon2)
            if (outmask & Caps.AZIMUTH):
                r.set_(azi12=_atan2d(lon12, psi12))
            if (outmask & Caps.DISTANCE):
                h = hypot(lon12, psi12)
                d = self._DIsometric2Rectifying(psi2, psi1)
                r.set_(s12=h * d * E.L / _90_0)
            if (outmask & Caps.AREA):
                r.set_(S12=self._S12d(lon12, psi2, psi1))
            if ((outmask | self._debug) & Caps._DEBUG_INVERSE):  # PYCHOK no cover
                r.set_(a=E.a, f=E.f, f1=E.f1, L=E.L,
                       b=E.b, e=E.e, e2=E.e2, k2=self._eF.k2,
                       lon12=lon12, psi1=psi1, exact=self.exact,
                       psi12=psi12, psi2=psi2)
        return r

    def Inverse7(self, lat1, lon1, azi12, s12, outmask=Caps.AZIMUTH_DISTANCE_AREA):
        '''Like method L{Rhumb.Inverse} but returning a L{Rhumb7Tuple} with area C{S12}.
        '''
        return self.Inverse(lat1, lon1, azi12, s12, outmask=outmask).toRhumb7Tuple()

    def InverseLine(self, lat1, lon1, lat2, lon2, name=NN):  # caps=Caps.STANDARD
        '''Define a L{RhumbLine} in terms of the I{inverse} rhumb problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).

           @return: A L{RhumbLine} instance and invoke its method
                    L{RhumbLine.Position} to compute each point.

           @note: Updates to this rhumb are reflected in the returned
                  rhumb line.
        '''
        r = self.Inverse(lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH)
        return RhumbLine(self, lat1=lat1, lon1=lon1, azi12=r.azi12,
                               caps=self._debug, name=name or self.name)

    Line = DirectLine  # synonimous

    def _MeanSinXid(self, psix, psiy):  # psix, psiy in degrees
        x, y = radians(psix), radians(psiy)
        t = _sincosSeries(False, _gd(x), _gd(y), *self._Rc2)
        return _Dlog(cosh(x), cosh(y)) * _Dcosh(x, y) + t * _Dgd(x, y)

    def orders(self, RAorder=None, TMorder=None):
        '''Get and set the I{RAorder} and/or I{TMorder}.

           @kwarg RAorder: I{Rhumb Area} order (C{int},
                           4, 5, 6, 7 or 8).
           @kwarg TMorder: I{Transverse Mercator} order (C{int},
                           4, 5, 6, 7 or 8).

           @return: L{RhumbOrder2Tuple}C{(RAorder, TMorder)} with the
                    previous C{RAorder} and C{TMorder} setting.
        '''
        t = RhumbOrder2Tuple(self.RAorder, self.TMorder)
        if RAorder not in (None, self._mRA):
            self.RAorder = RAorder
        if TMorder not in (None, self._mTM):
            self.TMorder = TMorder
        return t

    @Property
    def RAorder(self):
        '''Get the I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mRA

    @RAorder.setter  # PYCHOK setter!
    def RAorder(self, order):
        '''Get the I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        n = self._Xorder(_RCoeffs, RAorder=order)
        if n != self.RAorder:
            _update_all_rls(self)
            self._mRA = n

    @Property_RO
    def _Rc2(self):
        # for WGS84: (0, -0.0005583633519275459, -3.743803759172812e-07,  -4.633682270824446e-10,
        # RAorder 6:     -7.709197397676237e-13, -1.5323287106694307e-15, -3.462875359099873e-18)
        return self._Xs2(_RCoeffs, self.RAorder, RA=True)

    def _S12d(self, lon12, psi2, psi1):
        '''(INTERNAL) Compute the area C{S12}.
        '''
        r  = self._E.areax if self.exact else self._E.area
        # for WGS84: r = 510065621724088.44
        #            a = 6378137.0, f = 0.0033528106647474805
        #            b = 6356752.314245179, c2 = 40589732499314.76
        r *= lon12 / _720_0
        return r * self._MeanSinXid(psi2, psi1)

    @Property
    def TMorder(self):
        '''Get the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mTM

    @TMorder.setter  # PYCHOK setter!
    def TMorder(self, order):
        '''Set the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).

           @note: Setting C{TMorder} turns C{exact} off.
        '''
        n = self._Xorder(_AlpCoeffs, TMorder=order)
        if n != self.TMorder:
            _update_all_rls(self)
            self._mTM  = n
            self.exact = False

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
        return sep.join(pairs(d, prec=prec))

    def _Xorder(self, _Coeffs, **Xorder):
        '''(INTERNAL) Validate C{RAorder} or C{TMorder}.
        '''
        try:
            Xorder, order = Xorder.popitem()
            n = int(order)
            if n not in _Coeffs:
                t = sorted(map(str, _Coeffs.keys()))
                raise ValueError(_not_(_or(*t)))
        except (IndexError, TypeError, ValueError) as x:
            raise RhumbError(Xorder, order, txt=str(x))
        return n

    def _Xs2(self, _Coeffs, m, RA=False):
        '''(INTERNAL) Compute the C{R}, C{A} or C{B} terms for
           I{_sincosSeries}, return 2-tuple C{(terms, order)}.
        '''
        cs = _Coeffs[m]
        assert len(cs) == (((m + 1) * (m + 4)) if RA else
                            (m * (m + 3))) // 2
        n,  i = self._E.n, ((m + 2) if RA else 0)
        n_, X = n, [0]  # X[0] never used, it ...
        # ... is just an integration constant, so
        # it cancels when evaluating a definite
        # integral.  Don't bother computing it, it's
        # not used when invoking C{_sincosSeries}.
        for r in range(m - 1, -1, -1):  # [m-1 ... 0]
            j = i + r + 1
            X.append(_polynomial(n, cs, i, j) * n_ / cs[j])
            i   = j + 1
            n_ *= n
        return tuple(X), m


class RhumbError(_ValueError):
    '''Raised for an L{Rhumb} or L{RhumbLine} issue.
    '''
    pass


class _RhumbLine(_CapsBase):
    '''(INTERNAL) Class L{RhumbLine}
    '''
    _azi12     = _0_0
#   _caps      =  0
#   _lat1      = _0_0
#   _lon1      = _0_0
    _salp      = _0_0
    _calp      = _1_0
    _rhumb     =  None  # L{Rhumb} instance

    def __init__(self, rhumb, lat1, lon1, azi12, caps=0, name=NN):  # case=Caps.?
        '''New C{RhumbLine}.
        '''
        if (caps & Caps._DEBUG_DIRECT_LINE):
            self._debug |= caps & Caps._DEBUG_DIRECT_LINE
        _xinstanceof(Rhumb, rhumb=rhumb)
        self._lat1 = _Lat(lat1=_fix90(lat1))
        self._lon1 = _Lon(lon1=       lon1)
        if azi12:  # non-zero
            self.azi12 = azi12
        self._caps = caps
        if not (caps & Caps.LINE_OFF):
            _rls.append(self)
        if name:
            self.name=name or rhumb.name
        self._rhumb = rhumb  # last

    def __del__(self):  # XXX use weakref?
        if _rls:  # may be empty or None
            try:  # PYCHOK no cover
                _rls.remove(self)
            except (TypeError, ValueError):
                pass
        self._rhumb = None
        # _update_all(self)  # throws TypeError during Python 2 cleanup

    def _update(self, updated, *attrs, **unused):
        if updated:
            _update_all(self, *attrs)

    @Property
    def azi12(self):
        '''Get this rhumb line's I{azimuth} (compass C{degrees}).
        '''
        return self._azi12

    @azi12.setter  # PYCHOK setter!
    def azi12(self, azi12):
        '''Set this rhumb line's I{azimuth}.

           @arg azi12: The new I{azimuth} (compass C{degrees}).
        '''
        z = _norm180(azi12)
        if z != self._azi12:
            if self._rhumb:
                _update_all(self)
            self._azi12 = z
            self._salp, self._calp = sincos2d(z)  # no NEG0

    @Property_RO
    def ellipsoid(self):
        '''Get this rhumb line's ellipsoid (L{Ellipsoid}).
        '''
        return self._rhumb._E

    @property_RO
    def exact(self):
        '''Get this rhumb line's I{exact} option (C{bool}).
        '''
        return self.rhumb.exact

    @property_RO
    def lat1(self):
        '''Get this rhumb line's latitude (C{degrees90}).
        '''
        return self._lat1

    @property_RO
    def lon1(self):
        '''Get this rhumb line's longitude (C{degrees180}).
        '''
        return self._lon1

    @Property_RO
    def latlon1(self):
        '''Get this rhumb line's lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat1, self.lon1)

    @Property_RO
    def _mu1(self):
        '''(INTERNAL) Get the I{rectifying auxiliary} latitude C{mu} (C{degrees})..
        '''
        return self.ellipsoid.auxRectifying(self._lat1)

    @Property_RO
    def _psi1(self):
        '''(INTERNAL) Get the I{isometric auxiliary} latitude C{psi} (C{degrees}).
        '''
        return self.ellipsoid.auxIsometric(self._lat1)

    @property_RO
    def RAorder(self):
        '''Get this rhumb line's I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self.rhumb.RAorder

    @Property_RO
    def _r1rad(self):
        '''(INTERNAL) Get this rhumb line's I{icircle radius} (C{meter}).
        '''
        return radians(self.ellipsoid.circle4(self._lat1).radius)

    @Property_RO
    def rhumb(self):
        '''Get this rhumb line's rhumb (L{Rhumb}).
        '''
        return self._rhumb

    def Position(self, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Compute a position at a distance, optionally the area.

           @arg s12: The distance along this rhumb between its point and
                     the other point (C{meters}), can be negative.

           @return: L{GDict} with 2 up to 9 items C{azi12, s12, S12, lat2,
                    lon2, lat1, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the other point in C{degrees} and the area
                    C{S12} under the rhumb line in C{meter} I{squared}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses
                  a pole, the longitude of the second point is indeterminate
                  and C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{epsilon}**-2 (where C{epsilon} is 2**-52).
                  This position is extremely close to the actual pole and
                  allows the calculation to be carried out in finite terms.
        '''
        r = GDict(S12=NAN, name=self.name)
        if (outmask & Caps.LATITUDE_LONGITUDE_AREA):
            r.set_(s12=s12, azi12=self.azi12)
            E, R = self.ellipsoid, self.rhumb
            mu12 = s12  * self._calp * _90_0 / E.L
            mu2  = mu12 + self._mu1
            if fabs(mu2) > 90:  # PYCHOK no cover
                mu2 = _norm180(mu2)  # reduce to [-180, 180)
                if fabs(mu2) > 90:  # point on anti-meridian
                    mu2 = _norm180(_180_0 - mu2)
                lat2x = E.auxRectifying(mu2, inverse=True)
                lon2x = NAN
            else:
                psi2 = self._psi1
                if self._calp:
                    lat2x = E.auxRectifying(mu2, inverse=True)
                    psi12 = R._DRectifying2Isometric(radians(mu2),
                                                     radians(self._mu1)) * mu12
                    lon2x = psi12 * self._salp / self._calp
                    psi2 += psi12
                else:  # PYCHOK no cover
                    lat2x = self.lat1
                    lon2x = self._salp * s12 / self._r1rad

                if (outmask & Caps.AREA):
                    r.set_(S12=R._S12d(lon2x, self._psi1, psi2))

            if (outmask & Caps.LATITUDE):
                r.set_(lat2=lat2x, lat1=self.lat1)
            if (outmask & Caps.LONGITUDE):
                if (outmask & Caps.LONG_UNROLL) and not isnan(lat2x):
                    lon2x +=  self.lon1
                else:
                    lon2x  = _norm180(_norm180(self.lon1) + lon2x)
                r.set_(lon2=lon2x, lon1=self.lon1)
            if ((outmask | self._debug) & Caps._DEBUG_DIRECT_LINE):  # PYCHOK no cover
                r.set_(a=E.a, f=E.f, f1=E.f1, L=E.L, exact=R.exact,
                       b=E.b, e=E.e, e2=E.e2, k2=R._eF.k2,
                       calp=self._calp, mu1 =self._mu1,  mu12=mu12,
                       salp=self._salp, psi1=self._psi1, mu2=mu2)
        return r

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{RhumbLine} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: C{RhumbLine} (C{str}).
        '''
        d = dict(rhumb=self.rhumb, lat1=self.lat1, lon1=self.lon1,
                                   azi12=self.azi12, exact=self.exact)
        return sep.join(pairs(d, prec=prec))

    @property_RO
    def TMorder(self):
        '''Get this rhumb line's I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self.rhumb.TMorder


class RhumbLine(_RhumbLine):
    '''Compute one or more points on a single rhumb line.

       Class L{RhumbLine} facilitates the determination of points
       on a single rhumb line.  The starting point (C{lat1}, C{lon1})
       and the azimuth C{azi12} are specified once.  Calls to method
       L{RhumbLine.Position} return the location of an other point
       and optionally the distance C{s12} along the rhumb line and
       the corresponding area C{S12} under the rhumb line.
    '''
    def __init__(self, rhumb, lat1=0, lon1=0, azi12=None, caps=0, name=NN):  # case=Caps.?
        '''New L{RhumbLine}.

           @arg rhumb: The rhumb reference (L{Rhumb}).
           @kwarg lat1: The latitude of the starting point (C{degrees90}).
           @kwarg lon1: The longitude of the starting point (C{degrees180}).
           @kwarg azi12: The azimuth of this rhumb line (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities.  Use C{Caps.LINE_OFF} if updates
                        to the B{C{rhumb}} should I{not} be reflected in
                        this rhumb line.
           @kwarg name: Optional name (C{str}).
        '''
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            rhumb = rhumb.copy(deep=False, name=NN(_UNDER_, rhumb.name))
        _RhumbLine.__init__(self, rhumb, lat1, lon1, azi12, caps=caps, name=name)


class RhumbOrder2Tuple(_GTuple):
    '''2-Tuple C{(RAorder, TMorder)} with a I{Rhumb Area} and
       I{Transverse Mercator} order, both C{int}.
    '''
    _Names_ = (Rhumb.RAorder.name, Rhumb.TMorder.name)
    _Units_ = (      Int,                Int)


class Rhumb7Tuple(_GTuple):
    '''7-Tuple C{(lat1, lon1, lat2, lon2, azi12, s12, S12)} with lat- C{lat1}, C{lat2}
       and longitudes C{lon1_init__}, C{lon2} of both points, the azimuth of the rhumb line
       C{azi12}, the rhumb distance C{s12} and the area C{S12} under the rhumb line.
    '''
    _Names_ = (_lat1_, _lon1_, _lat2_, _lon2_, _azi12_, _s12_, _S12_)
    _Units_ = (_Lat,   _Lon,   _Lat,   _Lon,   _Azi,    _M,    _M2)

    def toDirect9Tuple(self, dflt=NAN, **azi1_azi2_m12_M12_M21):
        '''Convert this L{Rhumb7Tuple} result to a 9-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenDirect}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg azi1_azi2_m12_M12_M21: Optional keyword arguments
                  to specify or override L{Inverse10Tuple} items.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2, s12,
                    m12, M12, M21, S12)}
        '''
        d = dict(azi1=self.azi12, azi2=self.azi12, m12=self.s12, M12=_1_0, M21=_1_0)  # PYCHOK attr
        d.update(azi1_azi2_m12_M12_M21)
        return self._toTuple(Direct9Tuple, dflt, d)

    def toInverse10Tuple(self, dflt=NAN, **a12_m12_M12_M21_salp1_calp1_salp2_calp2):
        '''Convert this L{Rhumb7Tuple} to a 10-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenInverse}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg a12_m12_M12_M21_salp1_calp1_salp2_calp2: Optional keyword
                  arguments to specify or override L{Inverse10Tuple} items.

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1, salp2, calp2,
                    m12, M12, M21, S12)}.
        '''
        s, c = sincos2d(self.azi12)  # PYCHOK attr
        d = dict(a12=self.azi12, m12=self.s12, M12=_1_0, M21=_1_0,  # PYCHOK attr
                                 salp1=s, calp1=c, salp2=s, calp2=c)
        d.update(a12_m12_M12_M21_salp1_calp1_salp2_calp2)
        return self._toTuple(Inverse10Tuple, dflt, d)

    def _toTuple(self, nTuple, dflt, updates={}):
        '''(INTERNAL) Convert this C{Rhumb7Tuple} to an B{C{nTuple}}.
        '''
        r = self.toGDict(**updates)
        t = tuple(r.get(n, dflt) for n in nTuple._Names_)
        return nTuple(t, name=self.name)


# Use I{Divided Differences} to determine (mu2 - mu1) / (psi2 - psi1) accurately.
# Definition: _Df(x,y,d) = (f(x) - f(y)) / (x - y), @see W. M. Kahan & R. J.
# Fateman, "Symbolic computation of Divided Differences", SIGSAM Bull. 33(3),
# 7-28 (1999). U{ACM<https://DL.ACM.org/doi/pdf/10.1145/334714.334716> and
# U{UCB<https://www.CS.bBerkeley.edu/~fateman/papers/divdiff.pdf>}

def _Datan(x, y):
    xy = x * y
    r  = xy + _1_0
    d  = x - y
    if d:  # 2 * xy > -1 == 2 * xy + 1 > 0 == xy + r > 0 == xy > -r
        r = (atan(d / r) if xy > -r else (atan(x) - atan(y))) / d
    else:
        r = _1_0 / r
    return r


def _Dasinh(x, y):
    hx = hypot1(x)
    d  = x - y
    if d:
        hx *= y
        hy  = x * hypot1(y)
        t = (d * (x + y) / (hy + hx)) if (x * y) > 0 else (hy - hx)
        r = asinh(t) / d
    else:
        r = _1_0 / hx
    return r


def _Dcosh(x, y):
    return _Dsincos(x, y, sinh, sinh)


def _DeatanhE(x, y, E):
    # Deatanhe(x, y) = eatanhe((x - y) / (1 - e^2 * x * y)) / (x - y)
    e = _1_0 - E.e2 * x * y
    # assert not isnear0(e)
    d =  x - y
    return (E._es_atanh(d / e) / d) if d else (E.e2 / e)


def _DfEt(x, y, eF):  # x, y are tangents
    # eF = Elliptic(-E.e12)  # -E.e2 / (1 - E.e2)
    x, y, r = atan(x), atan(y), _1_0
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
        t2 =  t**2 + _1_0
        D *= _2_0 / t2
        s  =  D * d
        if s:
            c = -(t + _1_0) * (t - _1_0) / t2
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
    d = x - y
    # Changed atanh(t / (x + y)) to asinh(t / (2 * sqrt(x*y))) to avoid
    # taking atanh(1) when x is large and y is 1.  N.B. this routine is
    # invoked with positive x and y, so no need to guard against taking
    # the sqrt of a negative quantity.  This also fixes bogus results
    # being returned for the area when an endpoint is at a pole.
    return (asinh(d / (sqrt(x * y) * _2_0)) * _2_0 / d) if d else (_1_0 / x)


def _Dsin(x, y):
    return _Dsincos(x, y, sin, cos)


def _Dsincos(x, y, s_, c_):
    r = c_((x + y) * _0_5)
    d =    (x - y) * _0_5
    if d:
        r *= s_(d) / d
    return r


def _Dsinh(x, y):
    return _Dsincos(x, y, sinh, cosh)


def _Dtand(x, y):  # x, y in degrees
    tx, ty = _tand(x), _tand(y)
    txy = tx * ty
    r   = txy + _1_0
    d   = x - y
    if d:  # 2 * txy > -1 == 2 * txy + 1 > 0 == txy + r > 0 == txy > -r
        r = ((_tand(d) * r) if txy > -r else (tx - ty)) / radians(d)
    return r


def _gd(x):
    return atan(sinh(x))


def _sincosSeries(sinp, x, y, c, n):
    # N.B. n >= 0 and c[] has n+1 elements 0..n,
    #      of which c[0] is ignored.
    # Use Clenshaw summation to evaluate
    #   m = (g(x) + g(y)) / 2        -- mean value
    #   s = (g(x) - g(y)) / (x - y)  -- average slope
    # where
    #   g(x) = sum(c[j] * SC(2 * j * x), j = 1..n)
    #   SC = sinp ? sin : cos
    #   CS = sinp ? cos : sin
    # ...
    d = x - y
    sp, cp, sd, cd = sincos2_(x + y, d)
    m =  cp * cd * _2_0
    s = (sp * sd / d) if d else sp
    # 2x2 matrices in row-major order
    A0, A1  = m, (-s * d**2)
    A2, A3  = (-s * _4_0), m
    b2 = b1 = (_0_0,) * 4
    if n > 0:
        b1 = c[n], _0_0, _0_0, c[n]
    for j in range(n - 1, 0, -1):
        b1, b2, cj = b2, b1, c[j]  # c[0] unused
        # b1 = A * b2 - b1 + c[j] * I
        B0, B1, B2, B3 = b2
        b1 = (fsum1_(A0 * B0, A1 * B2, -b1[0], cj),
              fsum1_(A0 * B1, A1 * B3, -b1[1]),
              fsum1_(A2 * B0, A3 * B2, -b1[2]),
              fsum1_(A2 * B1, A3 * B3, -b1[3], cj))
    # Here are the full expressions for m and s
    # f01, f02, f11, f12 = (0, 0, cd * sp,  2 * sd * cp) if sinp else \
    #                      (1, 0, cd * cp, -2 * sd * sp)
    # m = -b2[1] * f02 + (c[0] - b2[0]) * f01 + b1[0] * f11 + b1[1] * f12
    # s = -b2[2] * f01 + (c[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12
    s = fsum1_(b1[2] * cd * sp,  b1[3] * sd * cp * _2_0) if sinp else \
        fsum1_(b1[2] * cd * cp, -b1[3] * sd * sp * _2_0, -b2[2])
    return s


# _Alp- and _BetCoeffs copied from I{Karney}'s U{TransverseMercator.cpp
# <https://GeographicLib.SourceForge.io/C++/doc/TransverseMercator_8cpp_source.html>}
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
_RCoeffs = {  # Generated by Maxima on 2015-05-15 08:24:04-04:00
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

assert set(_AlpCoeffs.keys()) == set(_BetCoeffs.keys())

__all__ += _ALL_DOCS(_RhumbLine)

if __name__ == '__main__':

    def _re(fmt, r3, x3):
        e3 = []
        for r, x in zip(r3, x3):
            e = abs(r - x) / abs(x)
            e3.append('%.g' % (e,))
        print((fmt % r3) + ' errors: ' + ', '.join(e3))

    # <https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>
    rhumb = Rhumb(exact=True)  # WGS84 default
    r = rhumb.Direct7(40.6, -73.8, 51, 5.5e6)  # from JFK about NE
    _re('# lat2=%.8f, lon2=%.8f, S12=%.1f', (r.lat2, r.lon2, r.S12), (71.68889988, 0.25551982, 44095641862956))
    r = rhumb.Inverse7(40.6, -73.8, 51.6, -0.5)  # JFK to LHR
    _re('# azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (77.76838971, 5771083.383, 37395209100030))
    r = rhumb.Inverse7(40.6, -73.8, 35.8, 140.3)  # JFK to Tokyo Narita
    _re('# azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (-92.38888798, 12782581.068, -63760642939073))

# % python3 -m pygeodesy.rhumbx

# lat2=71.68889988, lon2=0.25551982, S12=44054189889953.1 errors: 4e-11, 2e-08, 0.0009
# azi12=77.76838971, s12=5771083.383 S12=37362984299663.1 errors: 3e-12, 6e-11, 0.0009
# azi12=-92.38888798, s12=12782581.068 S12=-63665156875078.1 errors: 2e-11, 2e-11, 0.001

# **) MIT License
#
# Copyright (C) 2022-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
