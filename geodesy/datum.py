
# -*- coding: utf-8 -*-

# Python implementation of geodesy tools for an ellipsoidal earth model.
#
# Includes ellipsoid parameters and datums for different geographic coordinate
# systems and methods for converting between them and to cartesian coordinates.
# Transcribed from JavaScript originals by (C) Chris Veness 2005-2016 and
# published under the same MIT Licence, see <http://www.movable-type.co.uk/
# scripts/geodesy/docs/latlon-ellipsoidal.js.html>

# Historical geodetic datums: a latitude/longitude point defines a geographic
# location on or above/below the earth’s surface, measured in degrees from
# the equator and the International Reference Meridian and meters above the
# ellipsoid, and based on a given datum.  The datum is based on a reference
# ellipsoid and tied to geodetic survey reference points.

# Modern geodesy is generally based on the WGS84 datum (as used for instance
# by GPS systems), but previously various reference ellipsoids and datum
# references were used.

# This module extends the core latlon-ellipsoidal module to include ellipsoid
# parameters and datum transformation parameters, and methods for converting
# between different (generally historical) datums.

# The UK Ordnance Survey National Grid References are still based on the
# otherwise historical OSGB36 datum, q.v. Ordnance Survey 'A guide to
# coordinate systems in Great Britain' Section 6 <http://www.ordnancesurvey
# .co.uk/docs/support/guide-coordinate-systems-great-britain.pdf> and also
# <http://www.ordnancesurvey.co.uk/blog/2014/12/2>.

# This module is used for UK Ordnance Survey mapping, and for historical purposes.

from utils import fdot, fStr, radians

# all public contants, classes and functions
__all__ = ('R_KM', 'R_M', 'R_NM', 'R_SM',  # constants
           'Datum',  'Ellipsoid',  'Transform',  # classes
           'Datums', 'Ellipsoids', 'Transforms')  # enum-like
__version__ = '16.09.05'


class _Enum(dict):  # enum-like

    def __init__(self, name):
        self.name = name

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError("%s.%s doesn't exist" % (self.name, attr))

    def __repr__(self):
        return '\n'.join('%s.%s: %r' % (self.name, n, v) for n, v in sorted(self.items()))

    def __str__(self):
        return self.name + ', '.join(sorted('.' + n for n in self.keys()))

    def _assert(self, **kwds):
        # check names against given names
        for a, v in kwds.items():
            assert getattr(self, a) is v

Ellipsoids = _Enum('Ellipsoids')
Transforms = _Enum('Transforms')
Datums     = _Enum('Datums')


class _Base(object):

    name = ''

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self)

    def __str__(self):
        return self.toStr()  # PYCHOK expected

    def _fStr(self, prec, *attrs):
        t = fStr([getattr(self, a) for a in attrs], prec=prec, sep=' ')
        t = ['%s=%s' % (a, v) for a, v in zip(attrs, t.split())]
        return ', '.join(t + ['name=%r' % (self.name,)])

    def _register(self, enum, name):
        # add this enum name
        if name:
            if name in enum:
                raise NameError('%s.%s exists' % (enum.name, name))
            enum[name] = self
            self.name = name


class Ellipsoid(_Base):
    '''Ellipsoid with semi-major axis (a), semi-minor axis (b) and
       flattening (f) or inverse flattening (1/f).

       Flattening (f) is defined as a / (a − b), the 1st eccentricity
       squared (e2) as f * (2 - f) or (a^2 - b^2) / a^2, the 2nd
       eccentricity squared (e22) as e2 / (1 - e2) or (a^2 - b^2) / b^2
       and the mean radius (R) as (2 * a + b) / 3.
    '''
    a    = 0  # semi-major axis in meter
    b    = 0  # semi-minor axis in meter
    f    = 0  # inverse flattening: (a - b) / a
    # precomputed, frequently used values
    a2   = 0  # (1 / a^2) squared
    a2b2 = 1  # (a / b)^2 squared or 1 / (1 - f) squared
    e2   = 0  # 1st eccentricity squared: f * (2 - f)
    e4   = 0  # e2 squared
    e21  = 1  # (1 - e2)
    e22  = 0  # 2nd eccentricity squared: e2 / (1 - e2)
    R    = 0  # mean radius: (2 * a + b) / 3 per IUGG definition

    def __init__(self, a, b, f, name=''):
        self.a = float(a)  # major half-axis in meter
        self.b = float(b)  # minor half-axis in meter
        if f > 0:
            self.f = 1.0 / f   # inverse flattening
            self.e2 = self.f * (2 - self.f)  # 1st eccentricity squared
            self.e4 = self.e2 * self.e2  # for Nvector.Cartesian.toNvector
            self.e21 = 1 - self.e2  # for Nvector.Cartesian.toNvector
            self.e22 = self.e2 / self.e21  # 2nd eccentricity squared
            self.a2b2 = (self.a / self.b) ** 2  # for Nvector.toCartesian
        self.a2 = 1 / (self.a * self.a)  # for Nvector.Cartesian.toNvector
        self.R = (2 * self.a + self.b) / 3  # per IUGG definition for WGS84
        self._register(Ellipsoids, name)
        # some sanity checks to catch mistakes
        if self.a < self.b or min(self.a, self.b) < 1:
            raise AssertionError('%s=%0.9f vs %s=%0.9f' %
                                ('a', self.a, 'b', self.b))
        t = (self.a - self.b) / self.a
        if abs(self.f - t) > 1e-7:
            raise AssertionError('%s=%0.9f vs %s=%0.9f' %
                                ('1/f', self.f, '(a-b)/a', t))
        self._ab_90 = (self.a - self.b) / 90  # see radiusat below

    def __eq__(self, other):
        return self is other or (isinstance(other, Ellipsoid) and
                                 self.a == other.a and
                                 self.b == other.b)

    def radiusAt(self, lat):
        '''Approximate the ellipsoid radius at the given
           latitude in degrees by trivial interpolation.

           @param {degrees90} lat - Latitude in degrees.

           @returns {number} Radius at that latitude.
        '''
        # r = major - (major - minor) * |lat| / 90
        return self.a - self._ab_90 * min(90, abs(lat))

    def toStr(self, prec=8):
        '''Return this ellipsoid as a string.

           @param {int} [prec=8] - Number of decimals, unstripped.

           @returns {string} Ellipsoid string.
        '''
        return self._fStr(prec, 'a', 'b', 'f', 'e2', 'e22', 'R')


R_KM = 6371.008771415  # mean (spherical) earth radius in kilo meter
R_M  = R_KM * 1.0e3    # mean (spherical) earth radius in meter
R_NM = R_KM * 0.53996  # mean (spherical) earth radius in nautical miles
R_SM = R_KM * 0.62137  # mean (spherical) earth radius in statute miles
#R_? = 6372797.560856  # XXX some other mean radius?  # PYCHOK expected

Ellipsoids._assert(
    WGS84         = Ellipsoid(6378137.0,   6356752.31425, 298.257223563, 'WGS84'),
    Airy1830      = Ellipsoid(6377563.396, 6356256.909,   299.3249646,   'Airy1830'),
    AiryModified  = Ellipsoid(6377340.189, 6356034.448,   299.3249646,   'AiryModified'),
    Bessel1841    = Ellipsoid(6377397.155, 6356078.963,   299.152815351, 'Bessel1841'),
    Clarke1866    = Ellipsoid(6378206.4,   6356583.8,     294.978698214, 'Clarke1866'),
    Clarke1880IGN = Ellipsoid(6378249.2,   6356515.0,     293.466021294, 'Clarke1880IGN'),  # XXX confirm
    Intl1924      = Ellipsoid(6378388.0,   6356911.946,   297.0,         'Intl1924'),  # aka Hayford
    GRS80         = Ellipsoid(6378137.0,   6356752.31414, 298.257222101, 'GRS80'),
    WGS72         = Ellipsoid(6378135.0,   6356750.5,     298.26,        'WGS72'),
    Sphere        = Ellipsoid(R_M,         R_M,             0.0,         'Sphere'),  # pseudo
)


class Transform(_Base):
    '''Helmert transformation.
    '''
    tx = 0  # translate in meter
    ty = 0
    tz = 0

    rx = 0  # rotation in radians
    ry = 0
    rz = 0

    s  = 0  # ppm
    s1 = 1

    sx = 0  # rotation in degree seconds
    sy = 0
    sz = 0

    def __init__(self, name='', tx=0, ty=0, tz=0,
                                sx=0, sy=0, sz=0, s=0):
        if tx:
            self.tx = float(tx)
        if ty:
            self.ty = float(ty)
        if tz:
            self.tz = float(tz)
        if sx:  # secs to rads
            self.rx = radians(sx / 3600.0)
            self.sx = sx
        if sy:
            self.ry = radians(sy / 3600.0)
            self.sy = sy
        if sz:
            self.rz = radians(sz / 3600.0)
            self.sz = sz
        if s:
            self.s = float(s)
            self.s1 = s * 1.0e-6 + 1  # normalize ppm to (s + 1)
        self._register(Transforms, name)

    def __eq__(self, other):
        return self is other or (isinstance(other, Transform) and
                                 self.tx == other.tx and
                                 self.ty == other.ty and
                                 self.tz == other.tz and
                                 self.rx == other.rx and
                                 self.ry == other.ry and
                                 self.rz == other.rz and
                                 self.s  == other.s)

    def toStr(self, prec=4):
        '''Return this transform as a string.

           @param {int} [prec=4] - Number of decimals, unstripped.

           @returns {string} Transform string.
        '''
        return self._fStr(prec, 'tx', 'ty', 'tz',
                                'rx', 'ry', 'rz', 's', 's1',
                                'sx', 'sy', 'sz')

    def transform(self, x, y, z, inverse=False):
        '''Transform a (geocentric) Cartesian point, forward or inverse.
        '''
        if inverse:
            _xyz = -1, -x, -y, -z
        else:
            _xyz =  1,  x,  y,  z
        # x', y', z' = (.tx + x * .s1 - y * .rz + z * .ry,
        #               .ty + x * .rz + y * .s1 - z * .rx,
        #               .tz - x * .ry + y * .rx + z * .s1)
        return (fdot(_xyz, self.tx,  self.s1, -self.rz,  self.ry),
                fdot(_xyz, self.ty,  self.rz,  self.s1, -self.rx),
                fdot(_xyz, self.tz, -self.ry,  self.rx,  self.s1))


# <https://en.wikipedia.org/wiki/Helmert_transformation>
Transforms._assert(
    WGS84          = Transform('WGS84'),
    Clarke1866     = Transform('Clarke1866', tx=8, ty=-160, tz=-176),
    ED50           = Transform('ED50', tx=89.5, ty=93.8, tz=123.1,
                                                         sz=  0.156, s=-1.2),
    Irl1975        = Transform('Irl1975', tx=-482.530, ty=130.596, tz=-564.557,
                     # XXX rotation signs may be opposite, to be checked
                                          sx=  -1.042, sy= -0.214, sz=  -0.631,
                                           s=  -1.1),
    Krassovsky1940 = Transform('Krassovsky1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),
    MGI            = Transform('MGI', tx=-577.326, ty=-90.129, tz=-463.920,
                                      sx=   5.137, sy=  1.474, sz=   5.297,
                                       s=  -2.423),
    NAD27          = Transform('NAD27', tx=8, ty=-160, tz=-176),
    NAD83          = Transform('NAD83', tx= 1.004,  ty=-1.910,   tz=-0.515,
                                        sx= 0.0267, sy= 0.00034, sz= 0.011,
                                         s=-0.0015),
    NTF            = Transform('NTF', tx=-168, ty= -60, tz=320),  # XXX verify
    OSGB36         = Transform('OSGB36', tx=-446.448,  ty=125.157, tz=-542.060,
                                         sx=  -0.1502, sy=-0.2470, sz=-0.8421,
                                          s=  20.4894),
    TokyoJapan     = Transform('TokyoJapan', tx=148, ty=-507, tz=-685),
    WGS72          = Transform('WGS72', tz=-4.5, sz=0.554, s=-0.22),
)


class Datum(_Base):
    '''Create a new datum from the given ellipsoid and transform.
    '''
    def __init__(self, ellipsoid, transform=None, name=''):
        self.ellipsoid = ellipsoid or Ellipsoids.WGS84
        if not isinstance(self.ellipsoid, Ellipsoid):
            raise TypeError('%s not an %s: %r' % ('ellipsoid', Ellipsoid.__name__, self.ellipsoid))

        self.transform = transform or Transforms.WGS84
        if not isinstance(self.transform, Transform):
            raise TypeError('%s not a %s: %r' % ('transform', Transform.__name__, self.transform))

        self._register(Datums, name or self.transform.name or self.ellipsoid.name)

    def __eq__(self, other):
        return self is other or (isinstance(other, Datum) and
                                 self.ellipsoid == other.ellipsoid and
                                 self.transform == other.transform)

    def toStr(self, **unused):
        '''Return this datum as a string.

           @returns {string} Datum string.
        '''
        t = []
        for a in ('ellipsoid', 'transform'):
            v = getattr(self, a)
            t.append('%s=%ss.%s' % (a, v.__class__.__name__, v.name))
        return ', '.join(t)


# Datums with associated ellipsoid and Helmert transform parameters
# to convert from WGS84 into the given datum.  More are available at
# <http://earth-info.nga.mil/GandG/coordsys/datums/NATO_DT.pdf> and
# <http://www.fieldenmaps.info/cconv/web/cconv_params.js>.
Datums._assert(
    WGS84      = Datum(Ellipsoids.WGS84, Transforms.WGS84),

    # <http://www.icao.int/safety/pbn/documentation/eurocontrol/eurocontrol wgs 84 implementation manual.pdf>
    WGS72      = Datum(Ellipsoids.WGS72, Transforms.WGS72),

    # <http://www.geocachingtoolbox.com?page=datumEllipsoidDetails>
    TokyoJapan = Datum(Ellipsoids.Bessel1841, Transforms.TokyoJapan),

    # XXX psuedo-ellipsoids for spherical LatLon
    Sphere     = Datum(Ellipsoids.Sphere, Transforms.WGS84, name='Sphere'),

    # <http://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf>
    OSGB36     = Datum(Ellipsoids.Airy1830, Transforms.OSGB36),

    #  Nouvelle Triangulation Francaise (Paris)  XXX verify
    NTF        = Datum(Ellipsoids.Clarke1880IGN, Transforms.NTF),

    # NAD83 (2009) == WGS84 - <http://www.uvm.edu/giv/resources/WGS84_NAD83.pdf>
    # (If you *really* must convert WGS84<->NAD83, you need more than this!)
    NAD83      = Datum(Ellipsoids.GRS80, Transforms.NAD83),

    # <http://en.wikipedia.org/wiki/Helmert_transformation>
    NAD27      = Datum(Ellipsoids.Clarke1866, Transforms.NAD27),

    # <http://osi.ie/OSI/media/OSI/Content/Publications/transformations_booklet.pdf>
    Irl1975    = Datum(Ellipsoids.AiryModified, Transforms.Irl1975),

    # <http://www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
    ED50       = Datum(Ellipsoids.Intl1924, Transforms.ED50),
)


if __name__ == '__main__':

    from sys import argv
    from tests import Tests as _Tests

    class Tests(_Tests):

        def testDatum(self):
            E = Ellipsoid(1000, 1000, 0, name='TestEllipsiod')
            self.test('ellipsoid', E is Ellipsoids.TestEllipsiod, 'True')
#           print(Ellipsoid())

            T = Transform(name='TestTransform')
            self.test('transform', T is Transforms.TestTransform, 'True')
#           print(Transform())

            D = Datum(E, T, name='TestDatum')
            self.test('datum', D is Datums.TestDatum, 'True')
#           print(Datum())

            R, fmt = Ellipsoids.WGS84.R, '%.5f'
            self.test('meanR', R, fmt % (R_M,), fmt=fmt)

    t = Tests(__file__, __version__)
    t.testDatum()
    t.results()

    # print all
    if '-verbose' in argv or '-v' in argv:
        for e in (Ellipsoids, Transforms, Datums):
            print('\n%r' % (e,))

    # Typical test results (on MacOS X):

    # testing datum.py version 16.09.05
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all datum.py tests passed (Python 2.7.10)

# Ellipsoids.Airy1830: Ellipsoid(a=6377563.396, b=6356256.909, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370461.23366667, name='Airy1830')
# Ellipsoids.AiryModified: Ellipsoid(a=6377340.189, b=6356034.448, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370238.27533333, name='AiryModified')
# Ellipsoids.Bessel1841: Ellipsoid(a=6377397.155, b=6356078.963, f=0.00334277, e2=0.00667437, e22=0.00671922, R=6370291.091, name='Bessel1841')
# Ellipsoids.Clarke1866: Ellipsoid(a=6378206.4, b=6356583.8, f=0.00339008, e2=0.00676866, e22=0.00681478, R=6370998.86666667, name='Clarke1866')
# Ellipsoids.Clarke1880IGN: Ellipsoid(a=6378249.2, b=6356515.0, f=0.00340755, e2=0.00680349, e22=0.00685009, R=6371004.46666667, name='Clarke1880IGN')
# Ellipsoids.GRS80: Ellipsoid(a=6378137.0, b=6356752.31414, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77138, name='GRS80')
# Ellipsoids.Intl1924: Ellipsoid(a=6378388.0, b=6356911.946, f=0.003367, e2=0.00672267, e22=0.00676817, R=6371229.31533333, name='Intl1924')
# Ellipsoids.Sphere: Ellipsoid(a=6371008.771415, b=6371008.771415, f=0.0, e2=0.0, e22=0.0, R=6371008.771415, name='Sphere')
# Ellipsoids.TestEllipsiod: Ellipsoid(a=1000.0, b=1000.0, f=0.0, e2=0.0, e22=0.0, R=1000.0, name='TestEllipsiod')
# Ellipsoids.WGS72: Ellipsoid(a=6378135.0, b=6356750.5, f=0.00335278, e2=0.00669432, e22=0.00673943, R=6371006.83333333, name='WGS72')
# Ellipsoids.WGS84: Ellipsoid(a=6378137.0, b=6356752.31425, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77141667, name='WGS84')

# Transforms.Clarke1866: Transform(tx=8.0, ty=-160.0, tz=-176.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='Clarke1866')
# Transforms.ED50: Transform(tx=89.5, ty=93.8, tz=123.1, rx=0.0, ry=0.0, rz=0.0, s=-1.2, s1=1.0, sx=0.0, sy=0.0, sz=0.156, name='ED50')
# Transforms.Irl1975: Transform(tx=-482.53, ty=130.596, tz=-564.557, rx=-0.0, ry=-0.0, rz=-0.0, s=-1.1, s1=1.0, sx=-1.042, sy=-0.214, sz=-0.631, name='Irl1975')
# Transforms.Krassovsky1940: Transform(tx=-24.0, ty=123.0, tz=94.0, rx=-0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=-0.02, sy=0.26, sz=0.13, name='Krassovsky1940')
# Transforms.MGI: Transform(tx=-577.326, ty=-90.129, tz=-463.92, rx=0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=5.137, sy=1.474, sz=5.297, name='MGI')
# Transforms.NAD27: Transform(tx=8.0, ty=-160.0, tz=-176.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='NAD27')
# Transforms.NAD83: Transform(tx=1.004, ty=-1.91, tz=-0.515, rx=0.0, ry=0.0, rz=0.0, s=-0.0015, s1=1.0, sx=0.0267, sy=0.0003, sz=0.011, name='NAD83')
# Transforms.NTF: Transform(tx=-168.0, ty=-60.0, tz=320.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='NTF')
# Transforms.OSGB36: Transform(tx=-446.448, ty=125.157, tz=-542.06, rx=-0.0, ry=-0.0, rz=-0.0, s=20.4894, s1=1.0, sx=-0.1502, sy=-0.247, sz=-0.8421, name='OSGB36')
# Transforms.TestTransform: Transform(tx=0.0, ty=0.0, tz=0.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='TestTransform')
# Transforms.TokyoJapan: Transform(tx=148.0, ty=-507.0, tz=-685.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='TokyoJapan')
# Transforms.WGS72: Transform(tx=0.0, ty=0.0, tz=-4.5, rx=0.0, ry=0.0, rz=0.0, s=-0.22, s1=1.0, sx=0.0, sy=0.0, sz=0.554, name='WGS72')
# Transforms.WGS84: Transform(tx=0.0, ty=0.0, tz=0.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='WGS84')

# Datums.ED50: Datum(ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
# Datums.Irl1975: Datum(ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
# Datums.NAD27: Datum(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
# Datums.NAD83: Datum(ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
# Datums.NTF: Datum(ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
# Datums.OSGB36: Datum(ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
# Datums.Sphere: Datum(ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
# Datums.TestDatum: Datum(ellipsoid=Ellipsoids.TestEllipsiod, transform=Transforms.TestTransform)
# Datums.TokyoJapan: Datum(ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
# Datums.WGS72: Datum(ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
# Datums.WGS84: Datum(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)

    # testing datum.py version 16.09.05
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all datum.py tests passed (Python 3.5.1)

# Ellipsoids.Airy1830: Ellipsoid(a=6377563.396, b=6356256.909, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370461.23366667, name='Airy1830')
# Ellipsoids.AiryModified: Ellipsoid(a=6377340.189, b=6356034.448, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370238.27533333, name='AiryModified')
# Ellipsoids.Bessel1841: Ellipsoid(a=6377397.155, b=6356078.963, f=0.00334277, e2=0.00667437, e22=0.00671922, R=6370291.091, name='Bessel1841')
# Ellipsoids.Clarke1866: Ellipsoid(a=6378206.4, b=6356583.8, f=0.00339008, e2=0.00676866, e22=0.00681478, R=6370998.86666667, name='Clarke1866')
# Ellipsoids.Clarke1880IGN: Ellipsoid(a=6378249.2, b=6356515.0, f=0.00340755, e2=0.00680349, e22=0.00685009, R=6371004.46666667, name='Clarke1880IGN')
# Ellipsoids.GRS80: Ellipsoid(a=6378137.0, b=6356752.31414, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77138, name='GRS80')
# Ellipsoids.Intl1924: Ellipsoid(a=6378388.0, b=6356911.946, f=0.003367, e2=0.00672267, e22=0.00676817, R=6371229.31533333, name='Intl1924')
# Ellipsoids.Sphere: Ellipsoid(a=6371008.771415, b=6371008.771415, f=0.0, e2=0.0, e22=0.0, R=6371008.771415, name='Sphere')
# Ellipsoids.TestEllipsiod: Ellipsoid(a=1000.0, b=1000.0, f=0.0, e2=0.0, e22=0.0, R=1000.0, name='TestEllipsiod')
# Ellipsoids.WGS72: Ellipsoid(a=6378135.0, b=6356750.5, f=0.00335278, e2=0.00669432, e22=0.00673943, R=6371006.83333333, name='WGS72')
# Ellipsoids.WGS84: Ellipsoid(a=6378137.0, b=6356752.31425, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77141667, name='WGS84')

# Transforms.Clarke1866: Transform(tx=8.0, ty=-160.0, tz=-176.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='Clarke1866')
# Transforms.ED50: Transform(tx=89.5, ty=93.8, tz=123.1, rx=0.0, ry=0.0, rz=0.0, s=-1.2, s1=1.0, sx=0.0, sy=0.0, sz=0.156, name='ED50')
# Transforms.Irl1975: Transform(tx=-482.53, ty=130.596, tz=-564.557, rx=-0.0, ry=-0.0, rz=-0.0, s=-1.1, s1=1.0, sx=-1.042, sy=-0.214, sz=-0.631, name='Irl1975')
# Transforms.Krassovsky1940: Transform(tx=-24.0, ty=123.0, tz=94.0, rx=-0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=-0.02, sy=0.26, sz=0.13, name='Krassovsky1940')
# Transforms.MGI: Transform(tx=-577.326, ty=-90.129, tz=-463.92, rx=0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=5.137, sy=1.474, sz=5.297, name='MGI')
# Transforms.NAD27: Transform(tx=8.0, ty=-160.0, tz=-176.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='NAD27')
# Transforms.NAD83: Transform(tx=1.004, ty=-1.91, tz=-0.515, rx=0.0, ry=0.0, rz=0.0, s=-0.0015, s1=1.0, sx=0.0267, sy=0.0003, sz=0.011, name='NAD83')
# Transforms.NTF: Transform(tx=-168.0, ty=-60.0, tz=320.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='NTF')
# Transforms.OSGB36: Transform(tx=-446.448, ty=125.157, tz=-542.06, rx=-0.0, ry=-0.0, rz=-0.0, s=20.4894, s1=1.0, sx=-0.1502, sy=-0.247, sz=-0.8421, name='OSGB36')
# Transforms.TestTransform: Transform(tx=0.0, ty=0.0, tz=0.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='TestTransform')
# Transforms.TokyoJapan: Transform(tx=148.0, ty=-507.0, tz=-685.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='TokyoJapan')
# Transforms.WGS72: Transform(tx=0.0, ty=0.0, tz=-4.5, rx=0.0, ry=0.0, rz=0.0, s=-0.22, s1=1.0, sx=0.0, sy=0.0, sz=0.554, name='WGS72')
# Transforms.WGS84: Transform(tx=0.0, ty=0.0, tz=0.0, rx=0.0, ry=0.0, rz=0.0, s=0.0, s1=0.0, sx=0.0, sy=0.0, sz=0.0, name='WGS84')

# Datums.ED50: Datum(ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
# Datums.Irl1975: Datum(ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
# Datums.NAD27: Datum(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
# Datums.NAD83: Datum(ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
# Datums.NTF: Datum(ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
# Datums.OSGB36: Datum(ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
# Datums.Sphere: Datum(ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
# Datums.TestDatum: Datum(ellipsoid=Ellipsoids.TestEllipsiod, transform=Transforms.TestTransform)
# Datums.TokyoJapan: Datum(ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
# Datums.WGS72: Datum(ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
# Datums.WGS84: Datum(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
