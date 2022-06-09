
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ classes U{Rhumb
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>}.

Class L{RhumbLine} has been enhanced with methods C{intersection2} and C{nearestOn4} to find
the intersection of two rhumb lines, respectively the nearest point on a rumb line.

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

from pygeodesy.basics import copysign0, isnan, _xinstanceof
from pygeodesy.errors import IntersectionError, _ValueError, _xdatum, _xkwds
# from pygeodesy.etm import ExactTransverseMercator  # in ._eTM below
from pygeodesy.fmath import euclid, favg, hypot, hypot1
from pygeodesy.fsums import Fmt, fsum1_, pairs
from pygeodesy.interns import INT0, NAN, NN, PI_2, _azi12_, _COMMASPACE_, \
                             _convergence_, _ellipsoidal_, _EPSqrt as _TOL, \
                             _intersection_, _lat1_, _lat2_, _lon1_, _lon2_, \
                             _no_, _not_, _s12_, _S12_, _UNDER_, _0_0, _0_5, \
                             _1_0, _2_0, _4_0, _90_0, _180_0, _720_0  # PYCHOK used!
from pygeodesy.karney import _a12_, _atan2d, Caps, _CapsBase as _RhumbBase, \
                             _diff182, Direct9Tuple, _EWGS84, _fix90, GDict, \
                             _GTuple, Inverse10Tuple, _norm180, _tand
from pygeodesy.ktm import _Xellipsoid, _Xorder, _Xs, \
                          _AlpCoeffs, _BetCoeffs  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.namedTuples import Distance2Tuple, LatLon2Tuple, NearestOn4Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, property_RO, \
                           _update_all
# from pygeodesy.streprs import Fmt, pairs  # from .fsums
from pygeodesy.units import Bearing as _Azi, Degrees as _Deg, Int, Lat, Lon, \
                            Meter as _M, Meter2 as _M2
from pygeodesy.utily import sincos2_, sincos2d
# from pygeodesy.vector3d import intersection3d3, Vector3d  # in ._eTM3d, ... below

from math import asinh, atan, cos, cosh, fabs, radians, sin, sinh, sqrt

__all__ = _ALL_LAZY.rhumbx
__version__ = '22.06.08'

_rls   = []  # instances of C{RbumbLine} to be updated
_TRIPS = 65  # .intersection2, 18+


class _Lat(Lat):
    '''(INTERNAL) Latitude B{C{lat}}.
    '''
    def __init__(self, *lat, **Error_name):
        kwds = _xkwds(Error_name, clip=0, Error=RhumbError)
        Lat.__new__(_Lat, *lat, **kwds)


class _Lon(Lon):
    '''(INTERNAL) Longitude B{C{lon}}.
    '''
    def __init__(self, *lon, **Error_name):
        kwds = _xkwds(Error_name, clip=0, Error=RhumbError)
        Lon.__new__(_Lon, *lon, **kwds)


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


class Rhumb(_RhumbBase):
    '''Class to solve of the I{direct} and I{inverse rhumb} problems, accurately.

       @see: The U{Detailed Description<https://GeographicLib.SourceForge.io/C++/doc/
       classGeographicLib_1_1Rhumb.html>} of I{Karney}'s C++ C{Rhumb Class}.
    '''
    _E     = _EWGS84
    _exact =  True
    _mRA   =  6
    _mTM   =  6

    def __init__(self, a_earth=_EWGS84, f=None, exact=True, name=NN, **RA_TMorder):
        '''New L{Rhumb}.

           @kwarg a_earth: This rhumb's earth (L{Ellipsoid}, L{Ellipsoid2},
                           L{a_f2Tuple}, L{Datum}, 2-tuple C{(a, f)}) or the
                           (equatorial) radius (C{scalar}).
           @kwarg f: The ellipsoid's flattening (C{scalar}), iff B{C{a_earth}} is
                     a C{scalar}, ignored otherwise.
           @kwarg exact: If C{True}, use an addition theorem for elliptic integrals
                         to compute I{Divided differences}, otherwise use the Krüger
                         series expansion (C{bool}), see also property C{exact}.
           @kwarg name: Optional name (C{str}).
           @kwarg RA_TMorder: Optional keyword arguments B{C{RAorder}} and B{C{TMorder}}
                              to set the respective C{order}, see properties C{RAorder}
                              and C{TMorder} and method C{orders}.

           @raise RhumbError: Invalid B{C{a_earth}}, B{C{f}} or B{C{RA_TMorder}}.
        '''
        if f is not None:
            self.ellipsoid = a_earth, f
        elif a_earth not in (_EWGS84, None):
            self.ellipsoid = a_earth
        if not exact:
            self._exact = False
        if name:
            self.name = name
        if RA_TMorder:
            self.orders(**RA_TMorder)

    def _DConformal2Rectifying(self, chix, chiy):
        return _1_0 + _sincosSeries(True, chix, chiy,
                *self._DRectifying2ConformalA2)

    @Property_RO
    def _DRectifying2ConformalA2(self):
        m = self.TMorder
        return _Xs(_AlpCoeffs, m, self.ellipsoid), m

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

    @deprecated_method
    def Direct7(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):
        '''DEPRECATED, use method L{Rhumb.Direct8}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self.Direct8(lat1, lon1, azi12, s12, outmask=outmask)._to7Tuple()

    def Direct8(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):
        '''Like method L{Rhumb.Direct} but returning a L{Rhumb8Tuple} with area C{S12}.
        '''
        return self.Direct(lat1, lon1, azi12, s12, outmask=outmask).toRhumb8Tuple()

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

    def _DIsometric(self, latx, laty):  # degrees
        phix = radians(latx)
        phiy = radians(laty)
        return (_Dtand(latx, laty) * _Dasinh(_tand(latx), _tand(laty))  -
                 _Dsin(phix, phiy) * _DeatanhE(sin(phix), sin(phiy), self._E))

    def _DIsometric2Rectifying(self, psix, psiy):  # degrees
        if self.exact:
            E = self.ellipsoid
            latx = E.auxIsometric(psix, inverse=True)
            laty = E.auxIsometric(psiy, inverse=True)
            r = self._DRectifying(latx, laty) / self._DIsometric(latx, laty)
        else:
            x, y = radians(psix), radians(psiy)
            r = self._DConformal2Rectifying(_gd(x), _gd(y)) * _Dgd(x, y)
        return r

    def _DRectifying(self, latx, laty):  # degrees
        E, eF = self._E, self._eF
        tbetx = E.f1 * _tand(latx)
        tbety = E.f1 * _tand(laty)
        return (E.f1 * _DfEt(tbetx, tbety, eF) * E.b
                     * _Dtand(latx, laty) * PI_2
                     * _Datan(tbetx, tbety)) / E.L

    def _DRectifying2Conformal(self, mux, muy):  # radians
        return _1_0 - _sincosSeries(True, mux, muy,
                *self._DRectifying2ConformalB2)

    @Property_RO
    def _DRectifying2ConformalB2(self):
        m = self.TMorder
        return _Xs(_BetCoeffs, m, self.ellipsoid), m

    def _DRectifying2Isometric(self, mux, muy):  # degrees
        E = self.ellipsoid
        latx = E.auxRectifying(mux, inverse=True)
        laty = E.auxRectifying(muy, inverse=True)
        if self.exact:
            r = self._DIsometric(latx, laty) / self._DRectifying(latx, laty)
        else:
            x, y = radians(mux), radians(muy)
            r = self._DRectifying2Conformal(x, y) * \
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
           L{a_f2Tuple}, 2-tuple C{(a, f)}) or the (equatorial) radius (C{scalar}).
        '''
        E = _Xellipsoid(a_earth_f, RhumbError)
        if self._E != E:
            _update_all_rls(self)
            self._E = E

    @property_RO
    def equatoradius(self):
        '''Get the C{ellipsoid}'s equatorial radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    a = equatoradius

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
        if self._exact != x:
            _update_all_rls(self)
            self._exact = x

    def flattening(self):
        '''Get the C{ellipsoid}'s flattening (C{float}).
        '''
        return self.ellipsoid.f

    f = flattening

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH_DISTANCE):
        '''Solve the I{inverse rhumb} problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).

           @return: L{GDict} with 5 to 8 items C{lat1, lon2, lat2, lon2,
                    azi12, s12, S12, a12}, the azimuth C{azi12} of the rhumb
                    line in compass C{degrees} between C{-180} and C{+180},
                    the rhumb angle C{a12} and distance C{s12} between both
                    points in C{degrees} respectively C{meter} and the area
                    C{S12} under the rhumb line in C{meter} I{squared}.

           @note: The shortest rhumb line is found.  If the end points are
                  on opposite meridians, there are two shortest rhumb lines
                  and the East-going one is chosen.

                  If either point is a pole, the cosine of its latitude is
                  taken to be C{epsilon}**-2 (where C{epsilon} is 2.0**-52).
                  This position is extremely close to the actual pole and
                  allows the calculation to be carried out in finite terms.
        '''
        r = GDict(name=self.name)
        if (outmask & Caps.AZIMUTH_DISTANCE_AREA):
            r.set_(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2)
            E = self.ellipsoid
            psi1 = E.auxIsometric(lat1)
            psi2 = E.auxIsometric(lat2)
            psi12 = psi2 - psi1
            lon12, _ = _diff182(lon1, lon2)
            if (outmask & Caps.AZIMUTH):
                r.set_(azi12=_atan2d(lon12, psi12))
            if (outmask & Caps.DISTANCE):
                a12 = hypot(lon12, psi12) * self._DIsometric2Rectifying(psi2, psi1)
                s12 = a12 * E._L_90
                r.set_(s12=s12, a12=copysign0(a12, s12))
            if (outmask & Caps.AREA):
                r.set_(S12=self._S12d(lon12, psi2, psi1))
            if ((outmask | self._debug) & Caps._DEBUG_INVERSE):  # PYCHOK no cover
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

    def Inverse8(self, lat1, lon1, azi12, s12, outmask=Caps.AZIMUTH_DISTANCE_AREA):
        '''Like method L{Rhumb.Inverse} but returning a L{Rhumb8Tuple} with area C{S12}.
        '''
        return self.Inverse(lat1, lon1, azi12, s12, outmask=outmask).toRhumb8Tuple()

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

    Line = DirectLine  # synonymous

    def _MeanSinXid(self, psix, psiy):  # degrees
        x, y = radians(psix), radians(psiy)
        t = _sincosSeries(False, _gd(x), _gd(y), *self._RA2)
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
        '''Get the I{Rhumb Area} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        n = _Xorder(_RACoeffs, RhumbError, RAorder=order)
        if self._mRA != n:
            _update_all_rls(self)
            self._mRA = n

    def _S12d(self, lon12, psi2, psi1):  # degrees
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
        n = _Xorder(_AlpCoeffs, RhumbError, TMorder=order)
        if self._mTM != n:
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


class RhumbError(_ValueError):
    '''Raised for an L{Rhumb} or L{RhumbLine} issue.
    '''
    pass


class _RhumbLine(_RhumbBase):
    '''(INTERNAL) Class L{RhumbLine}
    '''
    _azi12 = _0_0
#   _lat1  = _0_0
#   _lon1  = _0_0
    _salp  = _0_0
    _calp  = _1_0
    _rhumb =  None  # L{Rhumb} instance

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

    def distance2(self, lat, lon):
        '''Return the distance and (initial) bearing of a point.

           @arg lat: Latitude of the point (C{degrees}).
           @arg lon: Longitude of the points (C{degrees}).

           @return: A L{Distance2Tuple}C{(distance, initial)} with the C{distance}
                    in C{meter} and C{initial} bearing in C{degrees}.

           @see: Methods L{RhumbLine.intersection2} and L{RhumbLine.nearestOn4}.
        '''
        r = self.rhumb.Inverse(self.lat1, self.lon1, lat, lon, Caps.AZIMUTH_DISTANCE)
        return Distance2Tuple(r.s12, r.azi12)

    @Property_RO
    def ellipsoid(self):
        '''Get this rhumb line's ellipsoid (L{Ellipsoid}).
        '''
        return self.rhumb._E

    @property_RO
    def exact(self):
        '''Get this rhumb line's I{exact} option (C{bool}).
        '''
        return self.rhumb.exact

    def intersection2(self, other, tol=_TOL, **eps):
        '''Find the intersection of this and an other rhumb line.

           @arg other: The other rhumb line ({RhumbLine}).
           @kwarg tol: Tolerance for longitudinal convergence (C{degrees}).
           @kwarg eps: Tolerance for L{intersection3d3} (C{EPS}).

           @return: A L{LatLon2Tuple}{(lat, lon)} with the C{lat}- and
                    C{lon}gitude of the intersection point.

           @raise IntersectionError: No convergence for this B{C{tol}} or
                                     no intersection for an other reason.

           @see: Methods L{RhumbLine.distance2} and L{RhumbLine.nearestOn4}
                 and function L{pygeodesy.intersection3d3}.

           @note: Each iteration involves a round trip to this rhumb line's
                  L{ExactTransverseMercator} or L{KTransverseMercator}
                  projection and invoking function L{intersection3d3} in
                  that domain.
        '''
        _xinstanceof(other, RhumbLine)
        _xdatum(self.rhumb, other.rhumb, Error=RhumbError)
        try:
            # set initial origin and halfway point
            p = LatLon2Tuple(favg(self.lat1, other.lat1),
                             favg(self.lon1, other.lon1))
            _xTM4 =  self.xTM.reverse  # ellipsoidal or spherical
            _i3d3 = _MODS.vector3d._intersect3d3  # NOT .intersection3d3!
            for i in range(1, _TRIPS):
                v = _i3d3(self._xTM3d(p),  self.azi12,
                         other._xTM3d(p), other.azi12, useZ=False, **eps)[0]
                t = _xTM4(v.x, v.y, lon0=p.lon)  # PYCHOK Reverse4Tuple
                d =  euclid(t.lon - p.lon, t.lat)  # PYCHOK t.lat + p.lat - p.lat
                p =  LatLon2Tuple(t.lat + p.lat, t.lon)  # PYCHOK p.
                if d < tol:
                    return LatLon2Tuple(p.lat, p.lon, iteration=i,  # PYCHOK p.
                                        name=self.intersection2.__name__)
        except Exception as x:
            raise IntersectionError(_no_(_intersection_), txt=str(x))
        t = _no_(_convergence_, Fmt.PAREN(d))
        raise IntersectionError(interation=i, tol=tol, txt=t)

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
        return self.ellipsoid.auxRectifying(self.lat1)

    def nearestOn4(self, lat, lon, tol=_TOL, **eps):
        '''Find the point on this rhumb line nearest to a given point.

           @arg lat: Latitude of the point (C{degrees}).
           @arg lon: Longitude of the point (C{degrees}).
           @kwarg tol: Longitudinal convergence tolerance (C{degrees}).
           @kwarg eps: Tolerance for L{intersection3d3} (C{EPS}).

           @return: A L{NearestOn4Tuple}C{(lat, lon, distance, normal)} with
                    the C{lat}- and C{lon}gitude of the nearest point on and
                    the C{distance} in C{meter} to this rhumb line and with the
                    azimuth of the C{normal}, perpendicular to this rhumb line.

           @raise IntersectionError: No convergence for this B{C{eps}} or
                                     no intersection for an other reason.

           @see: Methods L{RhumbLine.distance2} and L{RhumbLine.intersection2}
                 and function L{intersection3d3}.
        '''
        z = _norm180(self.azi12 + _90_0)  # perpendicular
        r = _RhumbLine(self.rhumb, lat, lon, z, caps=Caps.LINE_OFF)
        p =  self.intersection2(r, tol=tol, **eps)
        t =  r.distance2(p.lat, p.lon)
        return NearestOn4Tuple(p.lat, p.lon, t.distance, z,
                                             iteration=p.iteration)

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
        '''(INTERNAL) Get this rhumb line's I{icircle radius} (C{meter}).
        '''
        return radians(self.ellipsoid.circle4(self.lat1).radius)

    @Property_RO
    def rhumb(self):
        '''Get this rhumb line's rhumb (L{Rhumb}).
        '''
        return self._rhumb

    def Position(self, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Compute a position at a distance on this rhumb line.

           @arg s12: The distance along this rhumb between its point and
                     the other point (C{meters}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: L{GDict} with 4 to 8 items C{azi12, a12, s12, S12, lat2,
                    lon2, lat1, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the other point in C{degrees}, the rhumb angle
                    C{a12} between both points in C{degrees} and the area C{S12}
                    under the rhumb line in C{meter} I{squared}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses a
                  pole, the longitude of the second point is indeterminate and
                  C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{epsilon}**-2 (where C{epsilon} is 2**-52).
                  This position is extremely close to the actual pole and
                  allows the calculation to be carried out in finite terms.
        '''
        r = GDict(name=self.name)
        if (outmask & Caps.LATITUDE_LONGITUDE_AREA):
            E, R = self.ellipsoid, self.rhumb
            mu12 = s12  * self._calp / E._L_90
            mu2  = mu12 + self._mu1
            if fabs(mu2) > 90:  # PYCHOK no cover
                mu2 = _norm180(mu2)  # reduce to [-180, 180)
                if fabs(mu2) > 90:  # point on anti-meridian
                    mu2 = _norm180(_180_0 - mu2)
                lat2x = E.auxRectifying(mu2, inverse=True)
                lon2x = NAN
                if (outmask & Caps.AREA):
                    r.set_(S12=NAN)
            else:
                psi2 = self._psi1
                if self._calp:
                    lat2x = E.auxRectifying(mu2, inverse=True)
                    psi12 = R._DRectifying2Isometric(mu2,
                                               self._mu1) * mu12
                    lon2x = psi12 * self._salp / self._calp
                    psi2 += psi12
                else:  # PYCHOK no cover
                    lat2x =  self.lat1
                    lon2x =  self._salp * s12 / self._r1rad
#                   psi12 = _0_0
#               a12 = copysign0(hypot(psi12, lon2x - self.lon1), s12)
                if (outmask & Caps.AREA):
                    r.set_(S12=R._S12d(lon2x, self._psi1, psi2))
            r.set_(s12=s12, azi12=self.azi12, a12=s12 / E._L_90)
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
                                   azi12=self.azi12, exact=self.exact,
                                   TMorder=self.TMorder, xTM=self.xTM)
        return sep.join(pairs(d, prec=prec))

    @property_RO
    def TMorder(self):
        '''Get this rhumb line's I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self.rhumb.TMorder

    @Property_RO
    def xTM(self):
        '''Get this rhumb line's I{Transverse Mercator} projection (L{ExactTransverseMercator}
           if I{exact} and I{ellipsoidal}, L{KTransverseMercator} otherwise).
        '''
        E = self.ellipsoid
        return _MODS.etm.ExactTransverseMercator(E) if E.isEllipsoidal and self.exact else \
               _MODS.ktm.KTransverseMercator(E, TMorder=self.TMorder)

    def _xTM3d(self, latlon0, z=INT0):
        '''(INTERNAL) Convert this C{latlon1 - latlon0} to C{Vector3d}.
        '''
        t = self.xTM.forward(self.lat1 - latlon0.lat, self.lon1, lon0=latlon0.lon)
        return _MODS.vector3d.Vector3d(t.easting, t.northing, z)


class RhumbLine(_RhumbLine):
    '''Compute one or more points on a single rhumb line.

       Class L{RhumbLine} facilitates the determination of points on a
       single rhumb line.  The starting point (C{lat1}, C{lon1}) and the
       azimuth C{azi12} are specified once.  Method L{RhumbLine.Position}
       returns the location of an other point and optionally the distance
       C{s12} along the rhumb line and the corresponding area C{S12} under
       the rhumb line.

       The intersection of two rhumb lines can be found with method
       L{RhumbLine.intersection2}.  Similarly, method L{RhumbLine.nearestOn4}
       computes the nearest point on a rhumb line.
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


class Rhumb8Tuple(_GTuple):
    '''8-Tuple C{(lat1, lon1, lat2, lon2, azi12, s12, S12, a12)} with lat- C{lat1},
       C{lat2} and longitudes C{lon1}, C{lon2} of both points, the azimuth of the
       rhumb line C{azi12}, the distance C{s12}, the area C{S12} under the rhumb
       line and the angular distance C{a12} between both points.
    '''
    _Names_ = (_lat1_, _lon1_, _lat2_, _lon2_, _azi12_, _s12_, _S12_,  _a12_)
    _Units_ = (_Lat,   _Lon,   _Lat,   _Lon,   _Azi,    _M,    _M2,    _Deg)

    def toDirect9Tuple(self, dflt=NAN, **a12_azi1_azi2_m12_M12_M21):
        '''Convert this L{Rhumb8Tuple} result to a 9-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenDirect}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg a12_azi1_azi2_m12_M12_M21: Optional keyword arguments
                     to specify or override L{Inverse10Tuple} items.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2, s12,
                    m12, M12, M21, S12)}
        '''
        d = dict(azi1=self.azi12, azi2=self.azi12, m12=self.s12, M12=_1_0, M21=_1_0)  # PYCHOK attr
        d.update(a12_azi1_azi2_m12_M12_M21)
        return self._toTuple(Direct9Tuple, dflt, d)

    def toInverse10Tuple(self, dflt=NAN, **a12_m12_M12_M21_salp1_calp1_salp2_calp2):
        '''Convert this L{Rhumb8Tuple} to a 10-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenInverse}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg a12_m12_M12_M21_salp1_calp1_salp2_calp2: Optional keyword
                     arguments to specify or override L{Inverse10Tuple} items.

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1, salp2, calp2,
                    m12, M12, M21, S12)}.
        '''
        s, c = sincos2d(self.azi12)  # PYCHOK attr
        d = dict(m12=self.s12, M12=_1_0, M21=_1_0,  # PYCHOK attr
                 salp1=s, calp1=c, salp2=s, calp2=c)
        d.update(a12_m12_M12_M21_salp1_calp1_salp2_calp2)
        return self._toTuple(Inverse10Tuple, dflt, d)

    def _toTuple(self, nTuple, dflt, updates={}):
        '''(INTERNAL) Convert this C{Rhumb8Tuple} to an B{C{nTuple}}.
        '''
        r = self.toGDict(**updates)
        t = tuple(r.get(n, dflt) for n in nTuple._Names_)
        return nTuple(t, name=self.name)

    @deprecated_method
    def _to7Tuple(self):
        '''DEPRECATED, do not use!
        '''
        return _MODS.deprecated.Rhumb7Tuple(self[:-1])


# Use I{Divided Differences} to determine (mu2 - mu1) / (psi2 - psi1) accurately.
# Definition: _Df(x,y,d) = (f(x) - f(y)) / (x - y), @see W. M. Kahan & R. J.
# Fateman, "Symbolic computation of Divided Differences", SIGSAM Bull. 33(3),
# 7-28 (1999). U{ACM<https://DL.ACM.org/doi/pdf/10.1145/334714.334716>, @see
# U{UCB<https://www.CS.Berkeley.edu/~fateman/papers/divdiff.pdf>}, Dec 8, 1999.

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
        r =  asinh(t) / d
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
    d = x - y
    sp, cp, sd, cd = sincos2_(x + y, d)
    sd = (sd / d) if d else _1_0
    m = cp * cd * _2_0
    s = sp * sd
    # 2x2 matrices in row-major order
    A0, A1  = m, (-s * d**2)
    A2, A3  = (-s * _4_0), m
    b2 = b1 = (_0_0,) * 4
    if n > 0:
        b1 = C[n], _0_0, _0_0, C[n]
    for j in range(n - 1, 0, -1):
        b1, b2, Cj = b2, b1, C[j]  # C[0] unused
        # b1 = A * b2 - b1 + C[j] * I
        B0, B1, B2, B3 = b2
        b1 = (fsum1_(A0 * B0, A1 * B2, -b1[0], Cj),
              fsum1_(A0 * B1, A1 * B3, -b1[1]),
              fsum1_(A2 * B0, A3 * B2, -b1[2]),
              fsum1_(A2 * B1, A3 * B3, -b1[3], Cj))
    # Here are the full expressions for m and s
    # f01, f02, f11, f12 = (0, 0, cd * sp,  2 * sd * cp) if sinp else \
    #                      (1, 0, cd * cp, -2 * sd * sp)
    # m = -b2[1] * f02 + (C[0] - b2[0]) * f01 + b1[0] * f11 + b1[1] * f12
    # s = -b2[2] * f01 + (C[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12
    s = fsum1_(b1[2] * cd * sp,  b1[3] * sd * cp * _2_0) if sinp else \
        fsum1_(b1[2] * cd * cp, -b1[3] * sd * sp * _2_0, -b2[2])
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

__all__ += _ALL_DOCS(Caps, _RhumbLine)

if __name__ == '__main__':

    def _re(fmt, r3, x3):
        e3 = []
        for r, x in zip(r3, x3):
            e = abs(r - x) / abs(x)
            e3.append('%.g' % (e,))
        print((fmt % r3) + ' errors: ' + ', '.join(e3))

    # <https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>
    rhumb = Rhumb(exact=True)  # WGS84 default
    print('# %r\n' % rhumb)
    r = rhumb.Direct8(40.6, -73.8, 51, 5.5e6)  # from JFK about NE
    _re('# lat2=%.8f, lon2=%.8f, S12=%.1f', (r.lat2, r.lon2, r.S12), (71.68889988, 0.25551982, 44095641862956.148438))
    r = rhumb.Inverse8(40.6, -73.8, 51.6, -0.5)  # JFK to LHR
    _re('# azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (77.76838971, 5771083.383328, 37395209100030.367188))
    r = rhumb.Inverse8(40.6, -73.8, 35.8, 140.3)  # JFK to Tokyo Narita
    _re('# azi12=%.8f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), (-92.388887981699639, 12782581.0676841792, -63760642939072.492))

# % python3 -m pygeodesy.rhumbx

# Rhumb(RAorder=6, TMorder=6, ellipsoid=Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594), exact=True)

# lat2=71.68889988, lon2=0.25551982, S12=44095641862956.1 errors: 4e-11, 2e-08, 2e-15
# azi12=77.76838971, s12=5771083.383 S12=37395209100030.4 errors: 3e-12, 5e-15, 2e-16
# azi12=-92.38888798, s12=12782581.068 S12=-63760642939072.5 errors: 2e-16, 0, 0

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
