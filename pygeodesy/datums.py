
# -*- coding: utf-8 -*-

u'''Datums and transformations thereof.

Classes L{Datum} and L{Transform} and registries L{Datums} and
L{Transforms}, respectively.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
including datums and ellipsoid parameters for different geographic coordinate
systems and methods for converting between them and to cartesian coordinates.
Transcoded from JavaScript originals by I{(C) Chris Veness 2005-2016} and
published under the same MIT Licence**, see U{latlon-ellipsoidal.js
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

Historical geodetic datums: a latitude/longitude point defines a geographic
location on or above/below the earth’s surface, measured in degrees from
the equator and the International Reference Meridian and meters above the
ellipsoid, and based on a given datum.  The datum is based on a reference
ellipsoid and tied to geodetic survey reference points.

Modern geodesy is generally based on the WGS84 datum (as used for instance
by GPS systems), but previously various reference ellipsoids and datum
references were used.

The UK Ordnance Survey National Grid References are still based on the otherwise
historical OSGB36 datum, q.v. U{Ordnance Survey 'A guide to coordinate systems
in Great Britain', Section 6<https://www.OrdnanceSurvey.co.UK/docs/support/
guide-coordinate-systems-great-britain.pdf>}.

@var Datums.BD72: Datum(name='BD72', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.BD72)
@var Datums.DHDN: Datum(name='DHDN', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.DHDN)
@var Datums.ED50: Datum(name='ED50', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
@var Datums.GDA2020: Datum(name='GDA2020', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.GRS80: Datum(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.Irl1975: Datum(name='Irl1975', ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
@var Datums.Krassovski1940: Datum(name='Krassovski1940', ellipsoid=Ellipsoids.Krassovski1940, transform=Transforms.Krassovski1940)
@var Datums.Krassowsky1940: Datum(name='Krassowsky1940', ellipsoid=Ellipsoids.Krassowsky1940, transform=Transforms.Krassowsky1940)
@var Datums.MGI: Datum(name='MGI', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.MGI)
@var Datums.NAD27: Datum(name='NAD27', ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
@var Datums.NAD83: Datum(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
@var Datums.NTF: Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
@var Datums.OSGB36: Datum(name='OSGB36', ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
@var Datums.Potsdam: Datum(name='Potsdam', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.Bessel1841)
@var Datums.Sphere: Datum(name='Sphere', ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
@var Datums.TokyoJapan: Datum(name='TokyoJapan', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
@var Datums.WGS72: Datum(name='WGS72', ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
@var Datums.WGS84: Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)

@var Transforms.BD72: Transform(name='BD72', tx=106.86863, ty=-52.29778, tz=103.72389, rx=-0, ry=-0, rz=-0.00001, s=1.2727, s1=1, sx=-0.33657, sy=-0.45696, sz=-1.84218)
@var Transforms.Bessel1841: Transform(name='Bessel1841', tx=-582, ty=-105, tz=-414, rx=-0.00001, ry=-0, rz=0.00001, s=-8.3, s1=0.99999, sx=-1.04, sy=-0.35, sz=3.08)
@var Transforms.Clarke1866: Transform(name='Clarke1866', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.DHDN: Transform(name='DHDN', tx=-591.28, ty=-81.35, tz=-396.39, rx=0.00001, ry=-0, rz=-0.00001, s=-9.82, s1=0.99999, sx=1.477, sy=-0.0736, sz=-1.458)
@var Transforms.ED50: Transform(name='ED50', tx=89.5, ty=93.8, tz=123.1, rx=0, ry=0, rz=0, s=-1.2, s1=1, sx=0, sy=0, sz=0.156)
@var Transforms.Identity: Transform(name='Identity', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.Irl1965: Transform(name='Irl1965', tx=-482.53, ty=130.596, tz=-564.557, rx=0.00001, ry=0, rz=0, s=-8.15, s1=0.99999, sx=1.042, sy=0.214, sz=0.631)
@var Transforms.Irl1975: Transform(name='Irl1975', tx=-482.53, ty=130.596, tz=-564.557, rx=-0.00001, ry=-0, rz=-0, s=-1.1, s1=1, sx=-1.042, sy=-0.214, sz=-0.631)
@var Transforms.Krassovski1940: Transform(name='Krassovski1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.Krassowsky1940: Transform(name='Krassowsky1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.MGI: Transform(name='MGI', tx=-577.326, ty=-90.129, tz=-463.92, rx=0.00002, ry=0.00001, rz=0.00003, s=-2.423, s1=1, sx=5.137, sy=1.474, sz=5.297)
@var Transforms.NAD27: Transform(name='NAD27', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.NAD83: Transform(name='NAD83', tx=1.004, ty=-1.91, tz=-0.515, rx=0, ry=0, rz=0, s=-0.0015, s1=1, sx=0.0267, sy=0.00034, sz=0.011)
@var Transforms.NTF: Transform(name='NTF', tx=-168, ty=-60, tz=320, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.OSGB36: Transform(name='OSGB36', tx=-446.448, ty=125.157, tz=-542.06, rx=-0, ry=-0, rz=-0, s=20.4894, s1=1.00002, sx=-0.1502, sy=-0.247, sz=-0.8421)
@var Transforms.TokyoJapan: Transform(name='TokyoJapan', tx=148, ty=-507, tz=-685, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.WGS72: Transform(name='WGS72', tx=0, ty=0, tz=-4.5, rx=0, ry=0, rz=0, s=-0.22, s1=1, sx=0, sy=0, sz=0.554)
@var Transforms.WGS84: Transform(name='WGS84', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, neg_, _xinstanceof
from pygeodesy.ellipsoids import a_f2Tuple, _4Ecef, Ellipsoid, \
                                 Ellipsoid2, Ellipsoids, Vector3Tuple
# from pygeodesy.errors import _IsnotError  # from .fmath
from pygeodesy.fmath import fdot, fmean, Fmt, _IsnotError
from pygeodesy.interns import NN, _Airy1830_, _AiryModified_, _Bessel1841_, \
                             _cartesian_, _Clarke1866_, _Clarke1880IGN_, \
                             _COMMASPACE_, _DOT_, _ellipsoid_, _ellipsoidal_, \
                             _float as _F, _GRS80_, _Intl1924_, _Krassovski1940_, \
                             _NAD27_, _Krassowsky1940_, _NAD83_, _name_, _s_, \
                             _Sphere_, _spherical_, _sx_, _sy_, _sz_, _transform_, \
                             _tx_, _ty_, _tz_, _UNDER_, _WGS72_, _WGS84_, \
                             _0_0, _0_26, _1_0, _2_0, _8_0, _3600_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedEnum, _NamedEnumItem, \
                                    _lazyNamedEnumItem as _lazy, Property_RO
# from pygeodesy.namedTuples import Vector3Tuple  # from .ellipsoids
# from pygeodesy.props import Property_RO  # from .named
# from pygeodesy.streprs import Fmt  # from .fmath
from pygeodesy.units import Radius_

from math import radians

__all__ = _ALL_LAZY.datums
__version__ = '21.12.30'

_BD72_       = 'BD72'
_DHDN_       = 'DHDN'
_ED50_       = 'ED50'
_GDA2020_    = 'GDA2020'
_Identity_   = 'Identity'
_Inverse_    = 'Inverse'
_Irl1965_    = 'Irl1965'
_Irl1975_    = 'Irl1975'
_MGI_        = 'MGI'
_NTF_        = 'NTF'
_OSGB36_     = 'OSGB36'
_Potsdam_    = 'Potsdam'
_TokyoJapan_ = 'TokyoJapan'


def _r_s2(s_):
    '''(INTERNAL) rotation in C{radians} and C{degree seconds}.
    '''
    return _F(radians(s_ / _3600_0)), s_


class Transform(_NamedEnumItem):
    '''Helmert transformation.

       @see: L{Transform7Tuple}.
    '''
    tx = _0_0  # x translation (C{meter})
    ty = _0_0  # y translation (C{meter})
    tz = _0_0  # z translation (C{meter})

    rx = _0_0  # x rotation (C{radians})
    ry = _0_0  # y rotation (C{radians})
    rz = _0_0  # z rotation (C{radians})

    s  = _0_0  # scale ppm (C{float})
    s1 = _1_0  # scale + 1 (C{float})

    sx = _0_0  # x rotation (degree seconds)
    sy = _0_0  # y rotation (degree seconds)
    sz = _0_0  # z rotation (degree seconds)

    def __init__(self, name=NN, tx=0, ty=0, tz=0,
                                sx=0, sy=0, sz=0, s=0):
        '''New L{Transform}.

           @kwarg name: Optional, unique name (C{str}).
           @kwarg tx: Optional X translation (C{meter}).
           @kwarg ty: Optional Y translation (C{meter}).
           @kwarg tz: Optional Z translation (C{meter}).
           @kwarg s: Optional scale (C{float}), ppm.
           @kwarg sx: Optional X rotation (C{degree seconds}).
           @kwarg sy: Optional Y rotation (C{degree seconds}).
           @kwarg sz: Optional Z rotation (C{degree seconds}).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        if tx:
            self.tx = tx
        if ty:
            self.ty = ty
        if tz:
            self.tz = tz
        if sx:  # secs to rads
            self.rx, self.sx = _r_s2(sx)
        if sy:
            self.ry, self.sy = _r_s2(sy)
        if sz:
            self.rz, self.sz = _r_s2(sz)
        if s:
            self.s  =    s
            self.s1 = _F(s * 1e-6 + _1_0)  # normalize ppm to (s + 1)

        self._register(Transforms, name)

    def __eq__(self, other):
        '''Compare this and an other transform.

           @arg other: The other transform (L{Transform}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Transform) and
                                 self.tx == other.tx and
                                 self.ty == other.ty and
                                 self.tz == other.tz and
                                 self.rx == other.rx and
                                 self.ry == other.ry and
                                 self.rz == other.rz and
                                 self.s  == other.s)

    def __matmul__(self, other):  # PYCHOK Python 3.5+
        '''Helmert-transform cartesian B{C{other}}.

           @raise TypeError: Invalid B{C{other}}.
        '''
        try:  # only CartesianBase
            return other.toTransform(self)
        except AttributeError:
            pass
        raise _IsnotError(_cartesian_, other=other)

    def inverse(self, name=NN):
        '''Return the inverse of this transform.

           @kwarg name: Optional, unique name (C{str}).

           @return: Inverse (Transform).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        return Transform(name=name or (self.name + _Inverse_),
                         tx=-self.tx, ty=-self.ty, tz=-self.tz,
                         sx=-self.sx, sy=-self.sy, sz=-self.sz, s=-self.s)

    def toStr(self, prec=5):  # PYCHOK expected
        '''Return this transform as a string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Transform attributes (C{str}).
        '''
        return self._instr(prec, _tx_, _ty_, _tz_,
                                 'rx', 'ry', 'rz', _s_, 's1',
                                 _sx_, _sy_, _sz_)

    def transform(self, x, y, z, inverse=False):
        '''Transform a (geocentric) Cartesian point, forward or inverse.

           @arg x: X coordinate (C{meter}).
           @arg y: Y coordinate (C{meter}).
           @arg z: Z coordinate (C{meter}).
           @kwarg inverse: Optional direction, forward or inverse (C{bool}).

           @return: A L{Vector3Tuple}C{(x, y, z)}, transformed.
        '''
        if inverse:
            _xyz = neg_(_1_0, x, y, z)
            _s1  = self.s1 - _2_0  # == -(1 - s * 1e-6)) == -(1 - (s1 - 1))
        else:
            _xyz = _1_0, x, y, z
            _s1  =  self.s1
        # x', y', z' = (.tx + x * .s1 - y * .rz + z * .ry,
        #               .ty + x * .rz + y * .s1 - z * .rx,
        #               .tz - x * .ry + y * .rx + z * .s1)
        return Vector3Tuple(fdot(_xyz, self.tx,      _s1, -self.rz,  self.ry),
                            fdot(_xyz, self.ty,  self.rz,      _s1, -self.rx),
                            fdot(_xyz, self.tz, -self.ry,  self.rx,      _s1),
                            name=self.name)


class Transforms(_NamedEnum):
    '''(INTERNAL) L{Transform} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, **name_tx_ty_tz_sx_sy_sz_s):
        '''(INTERNAL) Instantiate the C{Transform}.
        '''
        return Transform(**name_tx_ty_tz_sx_sy_sz_s)

Transforms = Transforms(Transform)  # PYCHOK singleton
'''Some pre-defined L{Transform}s, all I{lazily} instantiated.'''
# <https://WikiPedia.org/wiki/Helmert_transformation> from WGS84
Transforms._assert(
    BD72           = _lazy(_BD72_, tx=_F(106.868628), ty=_F(-52.297783), tz=_F(103.723893),
                     # <https://www.NGI.Be/FR/FR4-4.shtm> ETRS89 == WG84
                     # <https://GeoRepository.com/transformation_15929/BD72-to-WGS-84-3.html>
                                    sx=_F(-0.33657),   sy=_F( -0.456955), sz=_F( -1.84218),
                                     s=_F( 1.2727)),
    Bessel1841     = _lazy(_Bessel1841_, tx=_F(-582.0),  ty=_F(-105.0),  tz=_F(-414.0),
                                         sx=_F(  -1.04), sy=_F(  -0.35), sz=_F(   3.08),
                                          s=_F(  -8.3)),
    Clarke1866     = _lazy(_Clarke1866_, tx=_F(8), ty=_F(-160), tz=_F(-176)),
    DHDN           = _lazy(_DHDN_, tx=_F(-591.28),  ty=_F(-81.35),   tz=_F(-396.39),
                                   sx=_F(   1.477), sy=_F( -0.0736), sz=_F(  -1.458),
                                    s=_F(  -9.82)),  # Germany
    ED50           = _lazy(_ED50_, tx=_F(89.5), ty=_F(93.8), tz=_F(123.1),
                     # <https://GeoNet.ESRI.com/thread/36583> sz=_F(-0.156)
                     # <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal.js>
                     # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
                                                             sz=_F(  0.156), s=_F(-1.2)),
    Identity       = _lazy(_Identity_),
    Irl1965        = _lazy(_Irl1965_, tx=_F(-482.530), ty=_F(130.596), tz=_F(-564.557),
                                      sx=_F(   1.042), sy=_F(  0.214), sz=_F(   0.631),
                                       s=_F(  -8.15)),
    Irl1975        = _lazy(_Irl1975_, tx=_F(-482.530), ty=_F(130.596), tz=_F(-564.557),
                     # XXX rotation signs may be opposite, to be checked
                                      sx=_F(  -1.042), sy=_F( -0.214), sz=_F(  -0.631),
                                       s=_F(  -1.1)),
    Krassovski1940 = _lazy(_Krassovski1940_, tx=_F(-24.0),  ty=_F(123.0), tz=_F(94.0),
                                             sx=_F( -0.02), sy=    _0_26, sz=_F( 0.13),
                                              s=_F( -2.423)),  # spelling
    Krassowsky1940 = _lazy(_Krassowsky1940_, tx=_F(-24.0),  ty=_F(123.0), tz=_F(94.0),
                                             sx=_F( -0.02), sy=    _0_26, sz=_F( 0.13),
                                              s=_F( -2.423)),  # spelling
    MGI            = _lazy(_MGI_, tx=_F(-577.326), ty=_F(-90.129), tz=_F(-463.920),
                                  sx=_F(   5.137), sy=_F(  1.474), sz=_F(   5.297),
                                   s=_F(  -2.423)),  # Austria
    NAD27          = _lazy(_NAD27_, tx=_8_0, ty=_F(-160), tz=_F(-176)),
    NAD83          = _lazy(_NAD83_, tx=_F( 1.004),  ty=_F(-1.910),   tz=_F(-0.515),
                                    sx=_F( 0.0267), sy=_F( 0.00034), sz=_F( 0.011),
                                     s=_F(-0.0015)),
    NTF            = _lazy(_NTF_, tx=_F(-168), ty=_F(-60), tz=_F(320)),  # XXX verify
    OSGB36         = _lazy(_OSGB36_, tx=_F(-446.448),  ty=_F(125.157),  tz=_F(-542.060),
                                     sx=_F(  -0.1502), sy=_F( -0.2470), sz=_F(  -0.8421),
                                      s=_F(  20.4894)),
    TokyoJapan     = _lazy(_TokyoJapan_, tx=_F(148), ty=_F(-507), tz=_F(-685)),
    WGS72          = _lazy(_WGS72_, tz=_F(-4.5), sz=_F(0.554), s=_F(-0.22)),
    WGS84          = _lazy(_WGS84_),  # unity
)


class Datum(_NamedEnumItem):
    '''Ellipsoid and transform parameters for an earth model.
    '''
    _ellipsoid = Ellipsoids.WGS84  # default ellipsoid (L{Ellipsoid}, L{Ellipsoid2})
    _transform = Transforms.WGS84  # default transform (L{Transform})

    def __init__(self, ellipsoid, transform=None, name=NN):
        '''New L{Datum}.

           @arg ellipsoid: The ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
           @kwarg transform: Optional transform (L{Transform}).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Datum with that B{C{name}} already exists.

           @raise TypeError: If B{C{ellipsoid}} is not an L{Ellipsoid}
                             nor L{Ellipsoid2} or B{C{transform}} is
                             not a L{Transform}.
        '''
        self._ellipsoid = ellipsoid or Datum._ellipsoid
        _xinstanceof(Ellipsoid, ellipsoid=self.ellipsoid)

        self._transform = transform or Datum._transform
        _xinstanceof(Transform, transform=self.transform)

        self._register(Datums, name or self.transform.name or self.ellipsoid.name)

    def __eq__(self, other):
        '''Compare this and an other datum.

           @arg other: The other datum (L{Datum}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Datum) and
                          self.ellipsoid == other.ellipsoid and
                          self.transform == other.transform)

    def __matmul__(self, other):  # PYCHOK Python 3.5+
        '''Convert cartesian or ellipsoidal B{C{other}} to this datum.

           @raise TypeError: Invalid B{C{other}}.
        '''
        try:  # only CartesianBase and EllipsoidalLatLonBase
            return other.toDatum(self)
        except AttributeError:
            pass
        raise _IsnotError(_cartesian_, _ellipsoidal_, other=other)

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{datum} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    @Property_RO
    def ellipsoid(self):
        '''Get this datum's ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return self._ellipsoid

    @Property_RO
    def exactTM(self):
        '''Get the C{ExactTM} projection (L{ExactTransverseMercator}).
        '''
        from pygeodesy.etm import ExactTransverseMercator
        return ExactTransverseMercator(datum=self)

    @Property_RO
    def isEllipsoidal(self):
        '''Check whether this datum is ellipsoidal (C{bool}).
        '''
        return self._ellipsoid.isEllipsoidal

    @Property_RO
    def isOblate(self):
        '''Check whether this datum's ellipsoidal is I{oblate} (C{bool}).
        '''
        return self._ellipsoid.isOblate

    @Property_RO
    def isProlate(self):
        '''Check whether this datum's ellipsoidal is I{prolate} (C{bool}).
        '''
        return self._ellipsoid.isProlate

    @Property_RO
    def isSpherical(self):
        '''Check whether this datum is (near-)spherical (C{bool}).
        '''
        return self._ellipsoid.isSpherical

    def toStr(self, **unused):  # PYCHOK expected
        '''Return this datum as a string.

           @return: Datum attributes (C{str}).
        '''
        t = [Fmt.EQUAL(_name_, repr(self.named))]
        for a in (_ellipsoid_, _transform_):
            v = getattr(self, a)
            t.append(NN(Fmt.EQUAL(a, v.classname), _s_, _DOT_, v.name))
        return _COMMASPACE_.join(t)

    @Property_RO
    def transform(self):
        '''Get this datum's transform (L{Transform}).
        '''
        return self._transform


def _En2(earth, name):
    '''(INTERNAL) Helper for C{_ellipsoid} and C{_ellipsoidal_datum}.
    '''
    if isinstance(earth, (Ellipsoid, Ellipsoid2)):
        E = earth
        n = NN(_UNDER_, name or E.name)
    elif isinstance(earth, Datum):
        E = earth.ellipsoid
        n = NN(_UNDER_, name or earth.name)
    elif isinstance(earth, a_f2Tuple):
        n = NN(_UNDER_, name or earth.name)
        E = Ellipsoid(earth.a, earth.b, name=n)
    elif isinstance(earth, (tuple, list)) and len(earth) == 2:  # no cover
        n = NN(_UNDER_, name or getattr(earth, _name_, NN))
        a_f = a_f2Tuple(earth)  # *earth
        E = Ellipsoid(a_f.a, a_f.b, name=n)  # PYCHOK .a
    else:
        E, n = None, NN
    return E, n


def _ellipsoid(earth, name=NN):  # in .trf
    '''(INTERNAL) Create an L{Ellipsoid} or L{Ellipsoid2} from L{Datum} or C{a_f2Tuple}.
    '''
    E, _ = _En2(earth, name)
    if not E:
        _xinstanceof(Ellipsoid, Ellipsoid2, a_f2Tuple, Datum, earth=earth)
    return E


def _ellipsoidal_datum(earth, name=NN, raiser=False):
    '''(INTERNAL) Create a L{Datum} from an L{Ellipsoid} or L{Ellipsoid2} or C{a_f2Tuple}.
    '''
    if isinstance(earth, Datum):
        d = earth
    else:
        E, n = _En2(earth, name)
        if not E:
            _xinstanceof(Datum, Ellipsoid, Ellipsoid2, a_f2Tuple, earth=earth)
        d = Datum(E, transform=Transforms.Identity, name=n)
    if raiser and not d.isEllipsoidal:
        raise _IsnotError(_ellipsoidal_, earth=earth)
    return d


def _mean_radius(radius, *lats):
    '''(INTERNAL) Compute the mean radius of a L{Datum} from an L{Ellipsoid},
       L{Ellipsoid2} or scalar earth C{radius} for several latitudes.
    '''
    if isscalar(radius):
        r =  Radius_(radius, low=0, Error=TypeError)
    else:
        E = _ellipsoidal_datum(radius).ellipsoid
        r =  fmean(map(E.Rgeocentric, lats)) if lats else E.Rmean
    return r


def _spherical_datum(earth, name=NN, raiser=False):
    '''(INTERNAL) Create a L{Datum} from an L{Ellipsoid}, L{Ellipsoid2} or scalar earth C{radius}.
    '''
    if isscalar(earth):
        n = NN(_UNDER_, name)
        r = Radius_(earth, Error=TypeError)
        E = Ellipsoid(r, r, name=n)
        d = Datum(E, transform=Transforms.Identity, name=n)
    else:
        d = _ellipsoidal_datum(earth, name=name)
    if raiser and not d.isSpherical:
        raise _IsnotError(_spherical_, earth=earth)
    return d


class Datums(_NamedEnum):
    '''(INTERNAL) L{Datum} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, ellipsoid_name, transform_name, name=NN):
        '''(INTERNAL) Instantiate the L{Datum}.
        '''
        return Datum(Ellipsoids.get(ellipsoid_name),
                     Transforms.get(transform_name), name=name)

Datums = Datums(Datum)  # PYCHOK singleton
'''Some pre-defined L{Datum}s, all I{lazily} instantiated.'''
# Datums with associated ellipsoid and Helmert transform parameters
# to convert from WGS84 into the given datum.  More are available at
# <https://Earth-Info.NGA.mil/GandG/coordsys/datums/NATO_DT.pdf> and
# <XXX://www.FieldenMaps.info/cconv/web/cconv_params.js>.
Datums._assert(
    # Belgian Datum 1972, based on Hayford ellipsoid.
    # <https://NL.WikiPedia.org/wiki/Belgian_Datum_1972>
    # <https://SpatialReference.org/ref/sr-org/7718/html/>
    BD72           = _lazy(_BD72_, _Intl1924_, _BD72_),

    # Germany <https://WikiPedia.org/wiki/Bessel-Ellipsoid>
    #         <https://WikiPedia.org/wiki/Helmert_transformation>
    DHDN           = _lazy(_DHDN_, _Bessel1841_, _DHDN_),

    # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
    ED50           = _lazy(_ED50_, _Intl1924_, _ED50_),

    # Australia <https://ICSM.Gov.AU/datum/gda2020-and-gda94-technical-manuals>
#   ADG66          = _lazy(_ADG66_, _ANS_, _WGS84_),  # XXX Transform?
#   ADG84          = _lazy(_ADG84_, _ANS_, _WGS84_),  # XXX Transform?
#   GDA94          = _lazy(_GDA94_, _GRS80_, _WGS84_),
    GDA2020        = _lazy(_GDA2020_, _GRS80_, _WGS84_),  # XXX Transform?

    # <https://WikiPedia.org/wiki/GRS_80>
    GRS80          = _lazy(_GRS80_, _GRS80_, _WGS84_),

    # <https://OSI.IE/OSI/media/OSI/Content/Publications/transformations_booklet.pdf>
#   Irl1975        = _lazy(_Irl1965_, _AiryModified_, _Irl1965_),
    Irl1975        = _lazy(_Irl1975_, _AiryModified_, _Irl1975_),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Krassovski1940 = _lazy(_Krassovski1940_, _Krassovski1940_, _Krassovski1940_),  # XXX spelling?
    Krassowsky1940 = _lazy(_Krassowsky1940_, _Krassowsky1940_, _Krassowsky1940_),  # XXX spelling?

    # Austria <https://DE.WikiPedia.org/wiki/Datum_Austria>
    MGI            = _lazy(_MGI_, _Bessel1841_, _MGI_),

    # <https://WikiPedia.org/wiki/Helmert_transformation>
    NAD27          = _lazy(_NAD27_, _Clarke1866_, _NAD27_),

    # NAD83 (2009) == WGS84 - <https://www.UVM.edu/giv/resources/WGS84_NAD83.pdf>
    # (If you *really* must convert WGS84<->NAD83, you need more than this!)
    NAD83          = _lazy(_NAD83_, _GRS80_, _NAD83_),

    #  Nouvelle Triangulation Francaise (Paris)  XXX verify
    NTF            = _lazy(_NTF_, _Clarke1880IGN_, _NTF_),

    # <https://www.OrdnanceSurvey.co.UK/docs/support/guide-coordinate-systems-great-britain.pdf>
    OSGB36         = _lazy(_OSGB36_, _Airy1830_, _OSGB36_),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Potsdam        = _lazy(_Potsdam_, _Bessel1841_, _Bessel1841_),

    # XXX psuedo-ellipsoids for spherical LatLon
    Sphere         = _lazy(_Sphere_, _Sphere_, _WGS84_),

    # <https://www.GeoCachingToolbox.com?page=datumEllipsoidDetails>
    TokyoJapan     = _lazy(_TokyoJapan_, _Bessel1841_, _TokyoJapan_),

    # <https://www.ICAO.int/safety/pbn/documentation/eurocontrol/eurocontrol%20wgs%2084%20implementation%20manual.pdf>
    WGS72          = _lazy(_WGS72_, _WGS72_, _WGS72_),

    WGS84          = _lazy(_WGS84_, _WGS84_, _WGS84_),
)

_WGS84 = Datums.WGS84  # PYCHOK exported internally

if __name__ == '__main__':

    from pygeodesy.interns import _COMMA_, _NL_, _NL_var_

    # __doc__ of this file, force all into registery
    for r in (Datums, Transforms):
        t = [NN] + r.toRepr(all=True).split(_NL_)
        print(_NL_var_.join(i.strip(_COMMA_) for i in t))

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
