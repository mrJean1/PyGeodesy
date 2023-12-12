
# -*- coding: utf-8 -*-

u'''(INTERNAL) base classes C{RhumbBase} and C{RhumbLineBase}, pure Python version of I{Karney}'s
C++ classes U{Rhumb<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>}
and U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>}
from I{GeographicLib versions 2.0} and I{2.2} and I{Karney}'s C++ example U{Rhumb intersect
<https://SourceForge.net/p/geographiclib/discussion/1026620/thread/2ddc295e/>}.

Class L{RhumbLineBase} has been enhanced with methods C{Intersecant2}, C{Intersection} and C{PlumbTo}
to iteratively find the intersection of a rhumb line and a circle or an other rhumb line, respectively
a perpendicular geodesic or other rhumb line.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, especially the U{Class List<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>},
the background information on U{Rhumb lines<https://GeographicLib.SourceForge.io/C++/doc/rhumb.html>},
the utily U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} and U{Online
rhumb line calculations<https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2014-2023) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, unsigned0, _xinstanceof
from pygeodesy.constants import EPS, EPS0, EPS1, INT0, NAN, _over, \
                               _EPSqrt as _TOL, _0_0, _0_01, _1_0, _90_0
from pygeodesy.datums import Datum, _earth_datum, _spherical_datum, _WGS84
from pygeodesy.errors import IntersectionError, itemsorted, RhumbError, \
                            _xdatum, _xkwds, _xkwds_pop, _Xorder
# from pygeodesy.etm import ExactTransverseMercator  # _MODS
from pygeodesy.fmath import euclid, favg, sqrt_a,  Fsum
# from pygeodesy.formy import opposing  # _MODS
# from pygeodesy.fsums import Fsum  # from .fmath
from pygeodesy.interns import NN, _coincident_, _COMMASPACE_, _Dash, \
                             _dunder_nameof, _parallel_, _too_, _under
from pygeodesy.karney import _atan2d, Caps, _CapsBase, _diff182, _fix90, \
                             _norm180, GDict
# from pygeodesy.ktm import KTransverseMercator, _AlpCoeffs  # _MODS
from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS
# from pygeodesy.named import notOverloaded  # _MODS
from pygeodesy.namedTuples import Distance2Tuple, LatLon2Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_RO, _update_all
from pygeodesy.streprs import Fmt, pairs
from pygeodesy.units import Float_, Lat, Lon, Meter, Radius_,  Int  # PYCHOK shared
from pygeodesy.utily import acos1, _azireversed, _loneg, sincos2d, sincos2d_, \
                           _unrollon, _Wrap
from pygeodesy.vector3d import _intersect3d3, Vector3d  # in .Intersection below

from math import cos, fabs

__all__ = ()
__version__ = '23.12.03'

_anti_ = _Dash('anti')
_rls   = []  # instances of C{RbumbLine...} to be updated
_TRIPS = 65  # .Intersection, .PlumbTo, 19+


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
       of any C{RhumbLine} instances tied to the given
       C{Rhumb} instance B{C{r}}.
    '''
    # _xinstanceof(_MODS.rhumb.aux_.RhumbAux, _MODS.rhumb.ekx.Rhumb, r=r)
    _update_all(r)
    for rl in _rls:  # PYCHOK use weakref?
        if rl._rhumb is r:
            _update_all(rl)


class RhumbBase(_CapsBase):
    '''(INTERNAL) Base class for C{rhumb.aux_.RhumbAux} and C{rhumb.ekx.Rhumb}.
    '''
    _datum = _WGS84
    _exact =  True
    _f_max = _0_01
    _mTM   =  6  # see .TMorder

    def __init__(self, a_earth, f, exact, name):
        '''New C{RhumbAux} or C{Rhumb}.
        '''
        _earth_datum(self, a_earth, f=f, name=name)
        if not exact:
            self.exact = False
        if name:
            self.name = name

    @Property_RO
    def a(self):
        '''Get the C{ellipsoid}'s equatorial radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    equatoradius = a

    def ArcDirect(self, lat1, lon1, azi12, a12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Solve the I{direct rhumb} problem, optionally with area.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees}).
           @arg a12: Angle along the rhumb line from the given to the
                     destination point (C{degrees}), can be negative.

           @return: L{GDict} with 2 up to 8 items C{lat2, lon2, a12, S12,
                    lat1, lon1, azi12, s12} with the destination point's
                    latitude C{lat2} and longitude C{lon2} in C{degrees},
                    the rhumb angle C{a12} in C{degrees} and area C{S12}
                    under the rhumb line in C{meter} I{squared}.

           @raise ImportError: Package C{numpy} not found or not installed,
                               only required for area C{S12} when C{B{exact}
                               is True} and L{RhumbAux}.

           @note: If B{C{a12}} is large enough that the rhumb line crosses
                  a pole, the longitude of the second point is indeterminate
                  and C{NAN} is returned for C{lon2} and area C{S12}.

           @note: If the given point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        s12 = a12 * self._mpd
        return self._DirectRhumb(lat1, lon1, azi12, a12, s12, outmask)

    @Property_RO
    def b(self):
        '''Get the C{ellipsoid}'s polar radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.b

    polaradius = b

    @property
    def datum(self):
        '''Get this rhumb's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this rhumb's datum (L{Datum}).

           @raise RhumbError: If C{abs(B{f}} exceeds non-zero C{f_max} and C{exact=False}.
        '''
        _xinstanceof(Datum, datum=datum)
        if self._datum != datum:
            self._exactest(self.exact, datum.ellipsoid, self.f_max)
            _update_all_rls(self)
            self._datum = datum

    def _Direct(self, ll1, azi12, s12, **outmask):
        '''(INTERNAL) Short-cut version, see .latlonBase.rhumb....
        '''
        return self.Direct(ll1.lat, ll1.lon, azi12, s12, **outmask)

    def Direct(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Solve the I{direct rhumb} problem, optionally with area.

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

           @raise ImportError: Package C{numpy} not found or not installed,
                               only required for area C{S12} when C{B{exact}
                               is True} and L{RhumbAux}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses
                  a pole, the longitude of the second point is indeterminate
                  and C{NAN} is returned for C{lon2} and area C{S12}.

           @note: If the given point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        a12 = _over(s12, self._mpd)
        return self._DirectRhumb(lat1, lon1, azi12, a12, s12, outmask)

    def Direct8(self, lat1, lon1, azi12, s12, outmask=Caps.LATITUDE_LONGITUDE_AREA):
        '''Like method L{Rhumb.Direct} but returning a L{Rhumb8Tuple} with area C{S12}.
        '''
        return self.Direct(lat1, lon1, azi12, s12, outmask=outmask).toRhumb8Tuple()

    def _DirectLine(self, ll1, azi12, **caps_name):
        '''(INTERNAL) Short-cut version, see .latlonBase.
        '''
        return self.DirectLine(ll1.lat, ll1.lon, azi12, **caps_name)

    def DirectLine(self, lat1, lon1, azi12, **caps_name):
        '''Define a C{RhumbLine} in terms of the I{direct} rhumb
           problem to compute several points on a single rhumb line.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees}).
           @kwarg caps_name: Optional keyword arguments C{B{name}=NN} and
                       C{B{caps}=Caps.STANDARD}, a bit-or'ed combination of
                       L{Caps} values specifying the required capabilities.
                       Include C{Caps.LINE_OFF} if updates to the B{C{rhumb}}
                       should I{not} be reflected in this rhumb line.

           @return: A C{RhumbLine...} instance and invoke its method
                    C{.Position} to compute each point.

           @note: Updates to this rhumb are reflected in the returned
                  rhumb line, unless C{B{caps} |= Caps.LINE_OFF}.
        '''
        return self._RhumbLine(self, lat1, lon1, azi12, **caps_name)

    Line = DirectLine  # synonyms

    def _DirectRhumb(self, lat1, lon1, azi12, a12, s12, outmask):
        '''(INTERNAL) See methods C{.ArcDirect} and C{.Direct}.
        '''
        rl = self._RhumbLine(self, lat1, lon1, azi12, caps=Caps.LINE_OFF,
                                                      name=self.name)
        return rl._Position(a12, s12, outmask | self._debug)  # lat2, lon2, S12

    @Property
    def ellipsoid(self):
        '''Get this rhumb's ellipsoid (L{Ellipsoid}).
        '''
        return self.datum.ellipsoid

    @ellipsoid.setter  # PYCHOK setter!
    def ellipsoid(self, a_earth_f):
        '''Set this rhumb's ellipsoid (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
           L{a_f2Tuple}) or (equatorial) radius and flattening (2-tuple C{(a, f)}).

           @raise RhumbError: If C{abs(B{f}} exceeds non-zero C{f_max} and C{exact=False}.
        '''
        self.datum = _spherical_datum(a_earth_f, Error=RhumbError)

    @Property
    def exact(self):
        '''Get the I{exact} option (C{bool}).
        '''
        return self._exact

    @exact.setter  # PYCHOK setter!
    def exact(self, exact):
        '''Set the I{exact} option (C{bool}).  If C{True}, use I{exact} rhumb
           expressions, otherwise a series expansion (accurate for oblate or
           prolate ellipsoids with C{abs(flattening)} below C{f_max}.

           @raise RhumbError: If C{B{exact}=False} and C{abs(flattening})
                              exceeds non-zero C{f_max}.

           @see: Option U{B{-s}<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}
                 and U{ACCURACY<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html#ACCURACY>}.
        '''
        x = bool(exact)
        if self._exact != x:
            self._exactest(x, self.ellipsoid, self.f_max)
            _update_all_rls(self)
            self._exact = x

    def _exactest(self, exact, ellipsoid, f_max):
        # Helper for property setters C{ellipsoid}, C{exact} and C{f_max}
        if fabs(ellipsoid.f) > f_max > 0 and not exact:
            raise RhumbError(exact=exact, f=ellipsoid.f, f_max=f_max)

    @Property_RO
    def f(self):
        '''Get the C{ellipsoid}'s flattening (C{float}).
        '''
        return self.ellipsoid.f

    flattening = f

    @property
    def f_max(self):
        '''Get the I{max.} flattening (C{float}).
        '''
        return self._f_max

    @f_max.setter  # PYCHOK setter!
    def f_max(self, f_max):  # PYCHOK no cover
        '''Set the I{max.} flattening, not to exceed (C{float}).

           @raise RhumbError: If C{exact=False} and C{abs(flattening})
                              exceeds non-zero C{f_max}.
        '''
        f = Float_(f_max=f_max, low=_0_0, high=EPS1)
        if self._f_max != f:
            self._exactest(self.exact, self.ellipsoid, f)
            self._f_max = f

    def _Inverse(self, ll1, ll2, wrap, **outmask):
        '''(INTERNAL) Short-cut version, see .latlonBase.rhumb....
        '''
        if wrap:
            ll2 = _unrollon(ll1, _Wrap.point(ll2))
        return self.Inverse(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **outmask)

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH_DISTANCE):
        '''Solve the I{inverse rhumb} problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).

           @return: L{GDict} with 4 to 9 items C{lat1, lon1, lat2, lon2,
                    azi12, azi21, s12, a12, S12}, the rhumb line's azimuth
                    C{azi12} and I{reverse} azimuth C{azi21}, both in
                    compass C{degrees} between C{-180} and C{+180}, the
                    rhumb distance C{s12} and rhumb angle C{a12} between
                    both points in C{meter} respectively C{degrees} and
                    the area C{S12} under the rhumb line in C{meter}
                    I{squared}.

           @raise ImportError: Package C{numpy} not found or not installed,
                               only required for L{RhumbAux} area C{S12}
                               when C{B{exact} is True}.

           @note: The shortest rhumb line is found.  If the end points are
                  on opposite meridians, there are two shortest rhumb lines
                  and the East-going one is chosen.

           @note: If either point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        r  = GDict(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2, name=self.name)
        Cs = Caps
        if (outmask & Cs.AZIMUTH_DISTANCE_AREA):
            lon12, _ = _diff182(lon1, lon2, K_2_0=True)
            y, x, s1, s2 = self._Inverse4(lon12, r, outmask)
            if (outmask & Cs.AZIMUTH):
                z = _atan2d(y, x)
                r.set_(azi12=z, azi21=_azireversed(z))
            if (outmask & Cs.AREA):
                S12 = self._S12d(s1, s2, lon12)
                r.set_(S12=unsigned0(S12))  # like .gx
        return r

    def _Inverse4(self, lon12, r, outmask):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, lon12, r, Caps.toStr(outmask))

    def Inverse8(self, lat1, lon1, azi12, s12, outmask=Caps.AZIMUTH_DISTANCE_AREA):
        '''Like method L{Rhumb.Inverse} but returning a L{Rhumb8Tuple} with area C{S12}.
        '''
        return self.Inverse(lat1, lon1, azi12, s12, outmask=outmask).toRhumb8Tuple()

    def _InverseLine(self, ll1, ll2, wrap, **caps_name):
        '''(INTERNAL) Short-cut version, see .latlonBase.
        '''
        if wrap:
            ll2 = _unrollon(ll1, _Wrap.point(ll2))
        return self.InverseLine(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **caps_name)

    def InverseLine(self, lat1, lon1, lat2, lon2, **caps_name):
        '''Define a C{RhumbLine} in terms of the I{inverse} rhumb problem.

           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg lat2: Latitude of the second point (C{degrees90}).
           @arg lon2: Longitude of the second point (C{degrees180}).
           @kwarg caps_name: Optional keyword arguments C{B{name}=NN} and
                       C{B{caps}=Caps.STANDARD}, a bit-or'ed combination of
                       L{Caps} values specifying the required capabilities.
                       Include C{Caps.LINE_OFF} if updates to the B{C{rhumb}}
                       should I{not} be reflected in this rhumb line.

           @return: A C{RhumbLine...} instance and invoke its method
                    C{ArcPosition} or C{Position} to compute points.

           @note: Updates to this rhumb are reflected in the returned
                  rhumb line, unless C{B{caps} |= Caps.LINE_OFF}.
        '''
        r = self.Inverse(lat1, lon1, lat2, lon2, outmask=Caps.AZIMUTH)
        return self._RhumbLine(self, lat1, lon1, r.azi12, **caps_name)

    @Property_RO
    def _mpd(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self)

    @property_RO
    def RAorder(self):
        '''Get the I{Rhumb Area} order, C{None} always.
        '''
        return None

    @property_RO
    def _RhumbLine(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, underOK=True)

    def _S12d(self, s1, s2, lon):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, s1, s2, lon)

    @Property
    def TMorder(self):
        '''Get the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self._mTM

    @TMorder.setter  # PYCHOK setter!
    def TMorder(self, order):
        '''Set the I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).

           @note: Setting C{TMorder} turns property C{exact} off, but only
                  for L{Rhumb} instances.
        '''
        m = _Xorder(_MODS.ktm._AlpCoeffs, RhumbError, TMorder=order)
        if self._mTM != m:
            _update_all_rls(self)
            self._mTM = m
            if self.exact and isinstance(self, _MODS.rhumb.ekx.Rhumb):
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
        return sep.join(pairs(itemsorted(d, asorted=False), prec=prec))


class RhumbLineBase(_CapsBase):
    '''(INTERNAL) Base class for C{rhumb.aux_.RhumbLineAux} and C{rhumb.ekx.RhumbLine}.
    '''
    _azi12 = _0_0
    _calp  = _1_0
#   _caps  =  0
#   _debug =  0
#   _lat1  = _0_0
#   _lon1  = _0_0
#   _lon12 = _0_0
    _Rhumb =  RhumbBase  # compatible C{Rhumb} class
    _rhumb =  None  # C{Rhumb} instance
    _salp  = _0_0
    _talp  = _0_0

    def __init__(self, rhumb, lat1, lon1, azi12, caps=Caps.STANDARD, name=NN):
        '''New C{RhumbLine} or C{RhumbLineAux}.
        '''
        _xinstanceof(self._Rhumb, rhumb=rhumb)

        self._lat1  = _Lat(lat1=_fix90(lat1))
        self._lon1  = _Lon(lon1=       lon1)
        self._lon12 = _norm180(self._lon1)
        if azi12:  # non-zero, non-None
            self.azi12 = _norm180(azi12)

        n = name or rhumb.name
        if n:
            self.name=n

        self._caps   =  caps
        self._debug |= (caps | rhumb._debug) & Caps._DEBUG_DIRECT_LINE
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            self._rhumb = rhumb.copy(deep=False, name=_under(rhumb.name))
        else:
            self._rhumb = rhumb
            _rls.append(self)

    def __del__(self):  # XXX use weakref?
        if _rls:  # may be empty or None
            try:  # PYCHOK no cover
                _rls.remove(self)
            except (TypeError, ValueError):
                pass
        self._rhumb = None
        # _update_all(self)  # throws TypeError during Python 2 cleanup

    def ArcPosition(self, a12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Compute a point at a given angular distance on this rhumb line.

           @arg a12: The angle along this rhumb line from its origin to the
                     point (C{degrees}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: L{GDict} with 4 to 8 items C{azi12, a12, s12, S12, lat2,
                    lon2, lat1, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the point in C{degrees}, the rhumb distance
                    C{s12} in C{meter} from the start point of and the area
                    C{S12} under this rhumb line in C{meter} I{squared}.

           @raise ImportError: Package C{numpy} not found or not installed,
                               only required for L{RhumbLineAux} area C{S12}
                               when C{B{exact} is True}.

           @note: If B{C{a12}} is large enough that the rhumb line crosses a
                  pole, the longitude of the second point is indeterminate and
                  C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        return self._Position(a12, self.degrees2m(a12), outmask)

    @Property
    def azi12(self):
        '''Get this rhumb line's I{azimuth} (compass C{degrees}).
        '''
        return self._azi12

    @azi12.setter  # PYCHOK setter!
    def azi12(self, azi12):
        '''Set this rhumb line's I{azimuth} (compass C{degrees}).
        '''
        z = _norm180(azi12)
        if self._azi12 != z:
            if self._rhumb:
                _update_all(self)
            self._azi12 = z
            self._salp, self._calp = t = sincos2d(z)  # no NEG0
            self._talp = _over(*t)

    @property_RO
    def azi12_sincos2(self):  # PYCHOK no cover
        '''Get the sine and cosine of this rhumb line's I{azimuth} (2-tuple C{(sin, cos)}).
        '''
        return self._scalp, self._calp

    @property_RO
    def datum(self):
        '''Get this rhumb line's datum (L{Datum}).
        '''
        return self.rhumb.datum

    def degrees2m(self, angle):
        '''Convert an angular distance along this rhumb line to C{meter}.

           @arg angle: Angular distance (C{degrees}).

           @return: Distance (C{meter}).
        '''
        return float(angle) * self.rhumb._mpd

    @deprecated_method
    def distance2(self, lat, lon):  # PYCHOK no cover
        '''DEPRECATED on 23.09.23, use method L{RhumbLineAux.Inverse} or L{RhumbLine.Inverse}.

           @return: A L{Distance2Tuple}C{(distance, initial)} with the C{distance}
                    in C{meter} and C{initial} bearing (azimuth) in C{degrees}.
        '''
        r = self.Inverse(lat, lon)
        return Distance2Tuple(r.s12, r.azi12)

    @property_RO
    def ellipsoid(self):
        '''Get this rhumb line's ellipsoid (L{Ellipsoid}).
        '''
        return self.rhumb.ellipsoid

    @property_RO
    def exact(self):
        '''Get this rhumb line's I{exact} option (C{bool}).
        '''
        return self.rhumb.exact

    def Intersecant2(self, lat0, lon0, radius, napier=True, **tol_eps):
        '''Compute the intersection(s) of this rhumb line and a circle.

           @arg lat0: Latitude of the circle center (C{degrees}).
           @arg lon0: Longitude of the circle center (C{degrees}).
           @arg radius: Radius of the circle (C{meter}, conventionally).
           @kwarg napier: If C{True}, apply I{Napier}'s spherical triangle
                          instead of planar trigonometry (C{bool}).
           @kwarg tol_eps: Optional keyword arguments, see method
                           method L{Intersection} for further details.

           @return: 2-Tuple C{(P, Q)} with both intersections (representing
                    a rhumb chord), each a L{GDict} from method L{Intersection}
                    extended to 18 items by C{lat3, lon3, azi03, a03, s03}
                    with azimuth C{azi03} of, distance C{a03} in C{degrees}
                    and C{s03} in C{meter} along the rhumb line from the circle
                    C{lat0, lon0} to the chord center C{lat3, lon3}.  If this
                    rhumb line is tangential to the circle, both points
                    are the same L{GDict} instance with distances C{s02} and
                    C{s03} near-equal to the B{C{radius}}.

           @raise IntersectionError: The circle and this rhumb line
                                     do not intersect.

           @raise UnitError: Invalid B{C{radius}}.
        '''
        r = Radius_(radius)
        p = q = self.PlumbTo(lat0, lon0, exact=None, **tol_eps)
        a = q.s02
        t = dict(lat3=q.lat2, lon3=q.lon2, azi03=q.azi02, a03=q.a02, s03=a)
        if a < r:
            t.update(iteration=q.iteration, lat0=q.lat1, lon0=q.lon1,  # or lat0, lon0
                     name=_dunder_nameof(self.Intersecant2, self.name))
            if fabs(a) < EPS0:  # coincident centers
                d, h = _0_0, r
            else:
                d = q.s12
                if napier:  # Napier rule (R1) cos(b) = cos(c) / cos(a)
                    # <https://WikiPedia.org/wiki/Spherical_trigonometry>
                    m =  self.rhumb._mpr
                    h = (acos1(cos(r / m) / cos(a / m)) * m) if m else _0_0
                else:
                    h = _copysign(sqrt_a(r, a), a)
            p = q = self.Position(d + h).set_(**t)
            if h:
                q = self.Position(d - h).set_(**t)
        elif a > r:
            t = _too_(Fmt.distant(a))
            raise IntersectionError(self, lat0, lon0, radius,
                                          txt=t, **tol_eps)
        else:  # tangential
            q.set_(**t)  # == p.set(_**t)
        return p, q

    def intersection2(self, other, **tol_eps):  # PYCHOK no cover
        '''DEPRECATED on 23.10.10, use method L{Intersection}.'''
        p = self.Intersection(other, **tol_eps)
        r = LatLon2Tuple(p.lat2, p.lon2, name=self.intersection2.__name__)
        r._iteration = p.iteration
        return r

    def Intersection(self, other, tol=_TOL, **eps):
        '''I{Iteratively} find the intersection of this and an other rhumb line.

           @arg other: The other rhumb line (C{RhumbLine}).
           @kwarg tol: Tolerance for longitudinal convergence and parallel
                       error (C{degrees}).
           @kwarg eps: Tolerance for L{pygeodesy.intersection3d3} (C{EPS}).

           @return: The intersection point, a L{Position}-like L{GDict} with
                    13 items C{lat1, lon1, azi12, a12, s12, lat2, lon2, lat0,
                    lon0, azi02, a02, s02, at} with the rhumb angle C{a02}
                    and rhumb distance C{s02} between the start point C{lat0,
                    lon0} of the B{C{other}} rhumb line and the intersection
                    C{lat2, lon2}, the azimuth C{azi02} of the B{C{other}}
                    rhumb line and the angle C{at} between both rhumb lines.
                    See method L{Position} for further details.

           @raise IntersectionError: No convergence for this B{C{tol}} or
                                     no intersection for an other reason.

           @see: Methods C{distance2} and C{PlumbTo} and function
                 L{pygeodesy.intersection3d3}.

           @note: Each iteration involves a round trip to this rhumb line's
                  L{ExactTransverseMercator} or L{KTransverseMercator}
                  projection and function L{pygeodesy.intersection3d3} in
                  that domain.
        '''
        _xinstanceof(RhumbLineBase, other=other)
        _xdatum(self.rhumb, other.rhumb, Error=RhumbError)
        try:
            if self.others(other) is self:
                raise ValueError(_coincident_)
            # make invariants and globals locals
            _s_3d, s_az =  self._xTM3d,  self.azi12
            _o_3d, o_az = other._xTM3d, other.azi12
            p = _MODS.formy.opposing(s_az, o_az, margin=tol)
            if p is not None:  # == t in (True, False)
                raise ValueError(_anti_(_parallel_) if p else _parallel_)
            _diff =  euclid  # approximate length
            _i3d3 = _intersect3d3  # NOT .vector3d.intersection3d3
            _LL2T =  LatLon2Tuple
            _xTMr =  self.xTM.reverse  # ellipsoidal or spherical
            # use halfway point as initial estimate
            p = _LL2T(favg(self.lat1, other.lat1),
                      favg(self.lon1, other.lon1))
            for i in range(1, _TRIPS):
                v = _i3d3(_s_3d(p), s_az,  # point + bearing
                          _o_3d(p), o_az, useZ=False, **eps)[0]
                t = _xTMr(v.x, v.y, lon0=p.lon)  # PYCHOK Reverse4Tuple
                d = _diff(t.lon - p.lon, t.lat)  # PYCHOK t.lat + p.lat - p.lat
                p = _LL2T(t.lat + p.lat, t.lon)  # PYCHOK t.lon + p.lon = lon0
                if d < tol:  # 19 trips
                    break
            else:
                raise ValueError(Fmt.no_convergence(d))

            n = _dunder_nameof(self.Intersection, self.name)
            r =  self.Inverse( p.lat, p.lon, outmask=Caps.DISTANCE)
            t =  other.Inverse(p.lat, p.lon, outmask=Caps.DISTANCE)
            P =  GDict(lat1=self.lat1, lat2=p.lat, lat0=other.lat1,
                       lon1=self.lon1, lon2=p.lon, lon0=other.lon1,
                       azi12= self.azi12, a12=r.a12, s12=r.s12,
                       azi02=other.azi12, a02=t.a12, s02=t.s12,
                       at=other.azi12 - self.azi12, name=n)
            P._iteration = i  # .set_(iteration=i, ...) only
        except Exception as x:
            raise IntersectionError(self, other, tol=tol,
                                          eps=eps, cause=x)
        return P

    def Inverse(self, lat2, lon2, wrap=False, **outmask):
        '''Return the rhumb angle, distance, azimuth, I{reverse} azimuth, etc. of
           a rhumb line between the given point and this rhumb line's start point.

           @arg lat2: Latitude of the point (C{degrees}).
           @arg lon2: Longitude of the points (C{degrees}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                        and B{C{lon2}} (C{bool}).

           @return: L{GDict} with 8 items C{a12, s12, azi12, azi21, lat1, lon1,
                    lat2, lon2}, the rhumb angle C{a12} and rhumb distance C{s12}
                    between both points in C{degrees} respectively C{meter}, the
                    rhumb line's azimuth C{azi12} and I{reverse} azimuth C{azi21}
                    both in compass C{degrees} between C{-180} and C{+180}.
        '''
        if wrap:
            _, lat2, lon2 = _Wrap.latlon3(self.lon1, _fix90(lat2), lon2, wrap)
        r = self.rhumb.Inverse(self.lat1, self.lon1,        lat2,  lon2, **outmask)
        return r

    @Property_RO
    def isLoxodrome(self):
        '''Is this rhumb line a meridional (C{None}), a parallel
           (C{False}) or a C{True} loxodrome?

           @see: I{Osborne's} U{2.5 Rumb lines and loxodromes
                 <https://Zenodo.org/record/35392>}, page 37.
        '''
        return bool(self._salp) if self._calp else None

    @Property_RO
    def lat1(self):
        '''Get this rhumb line's latitude (C{degrees90}).
        '''
        return self._lat1

    @Property_RO
    def lon1(self):
        '''Get this rhumb line's longitude (C{degrees180}).
        '''
        return self._lon1

    @Property_RO
    def latlon1(self):
        '''Get this rhumb line's lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat1, self.lon1)

    def m2degrees(self, distance):
        '''Convert a distance along this rhumb line to an angular distance.

           @arg distance: Distance (C{meter}).

           @return: Angular distance (C{degrees}).
        '''
        return _over(float(distance), self.rhumb._mpd)

    @property_RO
    def _mu1(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, underOK=True)

    def _mu2lat(self, mu2):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, mu2, underOK=True)

    @deprecated_method
    def nearestOn4(self, lat0, lon0, **exact_eps_est_tol):  # PYCHOK no cover
        '''DEPRECATED on 23.10.10, use method L{PlumbTo}.'''
        P =  self.PlumbTo(lat0, lon0, **exact_eps_est_tol)
        r = _MODS.deprecated.classes.NearestOn4Tuple(P.lat2, P.lon2, P.s12, P.azi02,
                                                     name=self.nearestOn4.__name__)
        r._iteration = P.iteration
        return r

    @deprecated_method
    def NearestOn(self, lat0, lon0, **exact_eps_est_tol):  # PYCHOK no cover
        '''DEPRECATED on 23.10.30, use method L{PlumbTo}.'''
        return self.PlumbTo(lat0, lon0, **exact_eps_est_tol)

    def PlumbTo(self, lat0, lon0, exact=None, eps=EPS, est=None, tol=_TOL):
        '''Compute the I{perpendicular} intersection of this rumb line with a geodesic
           from the given point, in part transcoded from I{Karney}'s C++ U{rhumb-intercept
           <https://SourceForge.net/p/geographiclib/discussion/1026620/thread/2ddc295e/>}.

           @arg lat0: Latitude of the point (C{degrees}).
           @arg lon0: Longitude of the point (C{degrees}).
           @kwarg exact: If C{None}, use a rhumb line perpendicular to this rhumb
                         line, otherwise use an I{exact} C{Geodesic...} from the
                         given point perpendicular to this rhumb line (C{bool} or
                         C{Geodesic...}), see method L{Ellipsoid.geodesic_}.
           @kwarg eps: Optional tolerance for L{pygeodesy.intersection3d3} (C{EPS}),
                       used only if C{B{exact} is None}.
           @kwarg est: Optional, initial estimate for the distance C{s12} of the
                       intersection I{along} this rhumb line (C{meter}), used only
                       if C{B{exact} is not None}.
           @kwarg tol: Longitudinal convergence tolerance (C{degrees}) or distance
                       tolerance (C(meter)) when C{B{exact} is None}, respectively
                       C{not None}.

           @return: The intersection point on this rhumb line, a L{GDict} from method
                    L{Intersection} if B{C{exact}=None}.  If B{C{exact}} is not C{None},
                    a L{Position}-like L{GDict} of 13 items C{azi12, a12, s12, lat2,
                    lat1, lat0, lon2, lon1, lon0, azi0, a02, s02, at} with distance
                    C{a02} in C{degrees} and C{s02} in C{meter} between the given point
                    C{lat0, lon0} and the intersection C{lat2, lon2}, geodesic azimuth
                    C{azi0} at the given point and the (perpendicular) angle C{at}
                    between the geodesic and this rhumb line at the intersection.  The
                    I{geodesic} azimuth at the intersection is C{(at + azi12)}.  See
                    method L{Position} for further details.

           @raise ImportError: I{Karney}'s U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               package not found or not installed.

           @raise IntersectionError: No convergence for this B{C{eps}} or no
                                     intersection for some other reason.

           @see: Methods C{distance2}, C{Intersecant2} and C{Intersection}
                 and function L{pygeodesy.intersection3d3}.
        '''
        Cs = Caps
        if exact is None:
            z  = _norm180(self.azi12 + _90_0)  # perpendicular azimuth
            rl =  RhumbLineBase(self.rhumb, lat0, lon0, z, caps=Cs.LINE_OFF)
            P  =  self.Intersection(rl, tol=tol, eps=eps)

        else:  # C{rhumb-intercept}
            E   = self.ellipsoid
            _gI = E.geodesic_(exact=exact).Inverse
            gm  = Cs.STANDARD | Cs._REDUCEDLENGTH_GEODESICSCALE  # ^ Cs.DISTANCE_IN
            if est is None:  # get an estimate from the "perpendicular" geodesic
                r = _gI(self.lat1, self.lon1, lat0, lon0, outmask=Cs.AZIMUTH_DISTANCE)
                d, _ = _diff182(r.azi2, self.azi12, K_2_0=True)
                _, s12 = sincos2d(d)
                s12 *= r.s12  # signed
            else:
                s12  = Meter(est=est)
            try:
                tol  = Float_(tol=tol, low=EPS, high=None)
                # def _over(p, q):  # see @note at method C{.Position}
                #     return (p / (q or _copysign(tol, q))) if isfinite(q) else NAN

                _ErT = E.rocPrimeVertical  # aka rocTransverse
                _S12 = Fsum(s12).fsum2_
                for i in range(1, _TRIPS):  # suffix 1 == C++ 2, 2 == C++ 3
                    P =  self.Position(s12)  # outmask = Cs.LATITUDE_LONGITUDE
                    r = _gI(lat0, lon0, P.lat2, P.lon2, outmask=gm)
                    d, _ = _diff182(self.azi12, r.azi2, K_2_0=True)
                    s, c, s2, c2 = sincos2d_(d, r.lat2)
                    c2 *= _ErT(r.lat2)
                    s  *= _over(s2 * self._salp, c2) - _over(s * r.M21, r.m12)
                    s12, t = _S12(c / s)  # XXX _over?
                    if fabs(t) < tol:  # or fabs(c) < EPS
                        break
                P.set_(azi0=r.azi1, a02=r.a12, s02=r.s12,  # azi2=r.azi2,
                       lat0=lat0, lon0=lon0, iteration=i, at=r.azi2 - self.azi12,
                       name=_dunder_nameof(self.PlumbTo, self.name))
            except Exception as x:  # Fsum(NAN) Value-, ZeroDivisionError
                raise IntersectionError(lat0, lon0, tol=tol, exact=exact,
                                        eps=eps, est=est, iteration=i, cause=x)

        return P

    def Position(self, s12, outmask=Caps.LATITUDE_LONGITUDE):
        '''Compute a point at a given distance on this rhumb line.

           @arg s12: The distance along this rhumb line from its origin to
                     the point (C{meters}), can be negative.
           @kwarg outmask: Bit-or'ed combination of L{Caps} values specifying
                           the quantities to be returned.

           @return: L{GDict} with 4 to 8 items C{azi12, a12, s12, S12, lat2,
                    lat1, lon2, lon1} with latitude C{lat2} and longitude
                    C{lon2} of the point in C{degrees}, the rhumb angle C{a12}
                    in C{degrees} from the start point of and the area C{S12}
                    under this rhumb line in C{meter} I{squared}.

           @raise ImportError: Package C{numpy} not found or not installed,
                               only required for L{RhumbLineAux} area C{S12}
                               when C{B{exact} is True}.

           @note: If B{C{s12}} is large enough that the rhumb line crosses a
                  pole, the longitude of the second point is indeterminate and
                  C{NAN} is returned for C{lon2} and area C{S12}.

                  If the first point is a pole, the cosine of its latitude is
                  taken to be C{sqrt(L{EPS})}.  This position is extremely
                  close to the actual pole and allows the calculation to be
                  carried out in finite terms.
        '''
        return self._Position(self.m2degrees(s12), s12, outmask)

    def _Position(self, a12, s12, outmask):
        '''(INTERNAL) C{Arc-/Position} helper.
        '''
        r  = GDict(azi12=self.azi12, a12=a12, s12=s12, name=self.name)
        Cs = Caps
        if (outmask & Cs.LATITUDE_LONGITUDE_AREA):
            if a12 or s12:
                mu12 = self._calp * a12
                mu2  = self._mu1 + mu12
                if fabs(mu2) > 90:  # past pole
                    mu2 = _norm180(mu2)  # reduce to [-180, 180)
                    if fabs(mu2) > 90:  # point on anti-meridian
                        mu2 = _norm180(_loneg(mu2))
                    lat2 = self._mu2lat(mu2)
                    lon2 = S12 = NAN
                else:
                    lat2, lon2, S1, S2 = self._Position4(a12, mu2, s12, mu12)
                    if (outmask & Cs.AREA):
                        S12 = self.rhumb._S12d(S1, S2, lon2)
                        S12 = unsigned0(S12)  # like .gx
#                   else:
#                       S12 = None  # unused
                    if (outmask & Cs.LONGITUDE):
                        if (outmask & Cs.LONG_UNROLL):
                            lon2 +=  self.lon1
                        else:
                            lon2  = _norm180(self._lon12 + lon2)
            else:  # coincident
                lat2, lon2 = self.latlon1
                S12 = _0_0

            if (outmask & Cs.AREA):
                r.set_(S12=S12)
            if (outmask & Cs.LATITUDE):
                r.set_(lat2=lat2, lat1=self.lat1)
            if (outmask & Cs.LONGITUDE):
                r.set_(lon2=lon2, lon1=self.lon1)
        return r

    def _Position4(self, a12, mu2, s12, mu12):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        _MODS.named.notOverloaded(self, a12, s12, mu2, mu12)

    @Property_RO
    def rhumb(self):
        '''Get this rhumb line's rhumb (L{RhumbAux} or L{Rhumb}).
        '''
        return self._rhumb

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
        return sep.join(pairs(itemsorted(d, asorted=False), prec=prec))

    @property_RO
    def TMorder(self):
        '''Get this rhumb line's I{Transverse Mercator} order (C{int}, 4, 5, 6, 7 or 8).
        '''
        return self.rhumb.TMorder

    @Property_RO
    def xTM(self):
        '''Get this rhumb line's I{Transverse Mercator} projection (L{ExactTransverseMercator}
           if I{exact} and I{ellipsoidal}, otherwise L{KTransverseMercator} for C{TMorder}).
        '''
        E = self.ellipsoid
        # ExactTransverseMercator doesn't handle spherical earth models
        return _MODS.etm.ExactTransverseMercator(E) if self.exact and E.isEllipsoidal else \
               _MODS.ktm.KTransverseMercator(E, TMorder=self.TMorder)

    def _xTM3d(self, latlon0, z=INT0, V3d=Vector3d):
        '''(INTERNAL) C{xTM.forward} this C{latlon1} to C{V3d} with B{C{latlon0}}
           as current intersection estimate and central meridian.
        '''
        t = self.xTM.forward(self.lat1 - latlon0.lat, self.lon1, lon0=latlon0.lon)
        return V3d(t.easting, t.northing, z)


class _PseudoRhumbLine(RhumbLineBase):
    '''(INTERNAL) Pseudo-rhumb line for a geodesic (line), see C{geodesicw._PlumbTo}.
    '''
    def __init__(self, gl, name=NN):
        R = RhumbBase(gl.geodesic.ellipsoid, None, True, name)
        RhumbLineBase.__init__(self, R, gl.lat1, gl.lon1, 0, caps=Caps.LINE_OFF)
        self._azi1 = self.azi12 = gl.azi1
        self._gl   = gl
        self._gD   = gl.geodesic.Direct

    def PlumbTo(self, lat0, lon0, **exact_eps_est_tol):  # PYCHOK signature
        P = RhumbLineBase.PlumbTo(self, lat0, lon0, **exact_eps_est_tol)
        P.set_(azi1=self._gl.azi1, azi2=_xkwds_pop(P, azi12=None))
        return P  # geodesic L{Position}

    def Position(self, s12, **unused):  # PYCHOK signature
        r = self._gD(self.lat1, self.lon1, self._azi1, s12)
        self._azi1     = r.azi1
        self.azi12 = z = r.azi2
        self._salp, _  = sincos2d(z)
        return r.set_(azi12=z)


__all__ += _ALL_DOCS(RhumbBase, RhumbLineBase)

if __name__ == '__main__':

    from pygeodesy import printf, Rhumb as Rh, RhumbAux as Ah
    from pygeodesy.basics import _zip
    from pygeodesy.ellipsoids import _EWGS84

    Al = Ah(_EWGS84).Line(30, 0, 45)
    Rl = Rh(_EWGS84).Line(30, 0, 45)

    for i in range(1, 10):
        s = .5e6 + 1e6 / i
        a = Al.Position(s).lon2
        r = Rl.Position(s).lon2
        e = (fabs(a - r) / a) if a else 0
        printf('# Position.lon2 %.14f vs %.14f, diff %g', r, a, e)

    for exact in (None, False, True):
        for est in (None, 1e6):
            a = Al.PlumbTo(60, 0, exact=exact, est=est)
            r = Rl.PlumbTo(60, 0, exact=exact, est=est)
            printf('# %s, iteration=%s, exact=%s, est=%s\n# %s, iteration=%s',
                   a.toRepr(), a.iteration, exact, est,
                   r.toRepr(), r.iteration, nl=1)

    NE_=(71.688899882813, 0.2555198244234, 44095641862956.11)
    LHR=(77.7683897102557, 5771083.38332803, 37395209100030.39)
    NRT=(-92.38888798169965, 12782581.067684170, -63760642939072.50)

    def _ref(fmt, r3, x3):
        e3 = []
        for r, x in _zip(r3, x3):  # strict=True
            e = fabs(r - x) / fabs(x)
            e3.append('%.g' % (e,))
        printf((fmt % r3) + ', rel errors: ' + ', '.join(e3))

    for R in (Ah, Rh):  # <https://GeographicLib.SourceForge.io/cgi-bin/RhumbSolve -p 9> version 2.2
        rh = R(exact=True)  # WGS84 default
        printf('# %r', rh, nl=1)
        r = rh.Direct8(40.6, -73.8, 51, 5.5e6)  # from JFK about NE
        _ref('# JFK NE lat2=%.12f, lon2=%.12f, S12=%.1f', (r.lat2, r.lon2, r.S12), NE_)
        r = rh.Inverse8(40.6, -73.8, 51.6, -0.5)  # JFK to LHR
        _ref('# JFK-LHR azi12=%.12f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), LHR)
        r = rh.Inverse8(40.6, -73.8, 35.8, 140.3)  # JFK to Tokyo Narita
        _ref('# JFK-NRT azi12=%.12f, s12=%.3f S12=%.1f', (r.azi12, r.s12, r.S12), NRT)

# % python3 -m pygeodesy.rhumb.bases

# Position.lon2 11.61455846901637 vs 11.61455846901637, diff 3.05885e-16
# Position.lon2 7.58982302826842 vs 7.58982302826842, diff 2.34045e-16
# Position.lon2 6.28526067416369 vs 6.28526067416369, diff 2.82623e-16
# Position.lon2 5.63938995325146 vs 5.63938995325146, diff 1.57495e-16
# Position.lon2 5.25385527435707 vs 5.25385527435707, diff 0
# Position.lon2 4.99764604290380 vs 4.99764604290380, diff 8.88597e-16
# Position.lon2 4.81503363740473 vs 4.81503363740473, diff 1.84459e-16
# Position.lon2 4.67828821748836 vs 4.67828821748835, diff 5.69553e-16
# Position.lon2 4.57205667906283 vs 4.57205667906283, diff 5.82787e-16

# Intersection(a02=17.798332, a12=19.521356, at=90.0, azi02=135.0, azi12=45.0, lat0=60.0, lat1=30.0, lat2=45.0, lon0=0.0, lon1=0.0, lon2=15.830286, name='Intersection', s02=1977981.142985, s12=2169465.957531), iteration=9, exact=None, est=None
# Intersection(a02=17.798332, a12=19.521356, at=90.0, azi02=135.0, azi12=45.0, lat0=60.0, lat1=30.0, lat2=45.0, lon0=0.0, lon1=0.0, lon2=15.830286, name='Intersection', s02=1977981.142985, s12=2169465.957531), iteration=9

# Intersection(a02=17.798332, a12=19.521356, at=90.0, azi02=135.0, azi12=45.0, lat0=60.0, lat1=30.0, lat2=45.0, lon0=0.0, lon1=0.0, lon2=15.830286, name='Intersection', s02=1977981.142985, s12=2169465.957531), iteration=9, exact=None, est=1000000.0
# Intersection(a02=17.798332, a12=19.521356, at=90.0, azi02=135.0, azi12=45.0, lat0=60.0, lat1=30.0, lat2=45.0, lon0=0.0, lon1=0.0, lon2=15.830286, name='Intersection', s02=1977981.142985, s12=2169465.957531), iteration=9

# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=5, exact=False, est=None
# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=5

# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=7, exact=False, est=1000000.0
# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=7

# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=5, exact=True, est=None
# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=5

# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=7, exact=True, est=1000000.0
# PlumbTo(a02=17.967658, a12=27.74256, at=90.0, azi0=113.73626, azi12=45.0, lat0=60, lat1=30.0, lat2=49.634582, lon0=0, lon1=0.0, lon2=25.767876, name='PlumbTo', s02=1997960.116871, s12=3083112.636236), iteration=7

# RhumbAux(RAorder=None, TMorder=6, ellipsoid=Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e21=0.99330562, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594), exact=True)
# JFK NE lat2=71.688899882813, lon2=0.255519824423, S12=44095641862956.1, rel errors: 4e-16, 2e-13, 4e-16
# JFK-LHR azi12=77.768389710256, s12=5771083.383 S12=37395209100030.4, rel errors: 5e-16, 3e-16, 8e-16
# JFK-NRT azi12=-92.388887981700, s12=12782581.068 S12=-63760642939072.5, rel errors: 0, 1e-16, 7e-16

# Rhumb(RAorder=6, TMorder=6, ellipsoid=Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e21=0.99330562, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594), exact=True)
# JFK NE lat2=71.688899882813, lon2=0.255519824423, S12=44095641862956.1, rel errors: 2e-16, 1e-13, 5e-16
# JFK-LHR azi12=77.768389710256, s12=5771083.383 S12=37395209100030.4, rel errors: 4e-16, 3e-16, 6e-16
# JFK-NRT azi12=-92.388887981700, s12=12782581.068 S12=-63760642939072.5, rel errors: 0, 1e-16, 1e-16

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
