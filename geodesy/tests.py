
# -*- coding: utf-8 -*-

# Tests for all geodesy modules.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# and <http://www.movable-type.co.uk/scripts/latlong.html>.

from datum import R_NM, Datums
from dms import F_D, F_DM, F_DMS, F_RAD, compassDMS, lonDMS, normDMS
from utils import degrees, fStr
from inspect import isclass, isfunction, ismethod, ismodule
from os.path import basename
import sys

__all__ = ('Tests',)
__version__ = '16.09.06'

try:
    _int = int, long
    _str = basestring
except NameError:  # Python 3+
    _int = int
    _str = str


def _type(obj, attr):
    t = getattr(obj, attr, None)
    if isclass(t):
        t = '() class'
    elif isinstance(t, float):
        t = ' float'
    elif isfunction(t):
        t = '() function'
    elif isinstance(t, _int):
        t = ' int'
    elif ismethod(t):
        t = '() method'
    elif ismodule(t):
        t = ' module'
    elif type(t) is property:
        t = ' property'
    elif isinstance(t, _str):
        t = ' str'
    else:
        t = ' attribute'
    return t


class Tests(object):
    '''Tests based on @examples from the original JavaScript code
       and examples in <http://williams.best.vwh.net/avform.htm>
    '''
    _file   = ''
    _name   = ''
    _prefix = '    # '
    _tests  = []

    failed = 0
    total  = 0

    def __init__(self, file, version):
        self._file = file
        self._name = basename(file)
        self.title(file, version)

    def printf(self, fmt, *args, **kwds):  # nl=0
        nl = '\n' * kwds.get('nl', 0)
        print(nl + self._prefix + (fmt % args))

    def results(self, nl=0):
        v = sys.version.split()[0]
        n = self.failed
        if n:
            p = '' if n == 1 else 's'
            r = '(%.1f%%) FAILED' % (100.0 * n / self.total)
        else:
            n, p, r = 'all', 's', 'passed'
        self.printf('%s %s test%s %s (Python %s)', n, self._name, p, r, v, nl=nl)

    def test(self, name, value, expect, fmt='%s'):
        self.total += 1  # tests
        f, v = '', fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect):
            self.failed += 1  # failures
            f = '  FAILED, expected %s' % (expect,)
        self.printf('test %d %s: %s%s', self.total, name, v, f)

    def title(self, module, version):
        self.printf('testing %s version %s', basename(module), version, nl=1)

    def testEllipsoidal(self, LatLon):
        p = LatLon(51.4778, -0.0016)
        d = p.convertDatum(Datums.OSGB36)
        self.test('convertDatum', d, '51.477284°N, 000.000020°E, -45.91m')  # 51.4773°N, 000.0000°E, -45.91m

    def testLatLon(self, LatLon):

        p = LatLon(52.20472, 0.14056)
        self.test('lat/lonDMS', p, '52.20472°N, 000.14056°E')  # 52.20472°N, 000.14056°E
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 3),  '''52°12.283'N, 000°08.434'E''')
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 4),  '''52°12.2832'N, 000°08.4336'E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 0), '''52°12'17"N, 000°08'26"E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 1), '''52°12'17.0"N, 000°08'26.0"E''')
        self.test('lat/lonDMS F_RAD', p.toStr(F_RAD, 6), '0.911144N, 0.002453E')
        q = LatLon(*map(degrees, p.toradians()))
        self.test('equals', q.equals(p), 'True')

        LAX = LatLon(33.+57./60, -(118.+24./60))
        JFK = LatLon(degrees(0.709186), -degrees(1.287762))

        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        self.test('equals', p.equals(q), 'False')

        b = p.bearingTo(q)
        self.test('bearingTo', b, '156.1666', '%.4f')  # 156.2
        b = p.finalBearingTo(q)
        self.test('finalBearingTo', b, '157.8904', '%.4f')
        b = LAX.bearingTo(JFK)
        self.test('bearingTo', b, '65.8921', '%.4f')  # 66

        c = p.copy()
        self.test('copy', p.equals(c), 'True')

        d = p.distanceTo(q)
        self.test('distanceTo', d, '404279.720589', '%.6f')  # 404300
        d = q.distanceTo(p)
        self.test('distanceTo', d, '404279.720589', '%.6f')  # 404300
        d = LAX.distanceTo(JFK, radius=R_NM)
        self.test('distanceTo', d, '2145', '%.0f')  # 2144

        m = p.midpointTo(q)
        self.test('midpointTo', m, '50.536327°N, 001.274614°E')  # 50.5363°N, 001.2746°E

        p = LatLon(51.4778, -0.0015)
        d = p.destination(7794, 300.7)
        self.test('destination', d, '51.513546°N, 000.098345°W')  # 51.5135°N, 0.0983°W ???
        self.test('destination', d.toStr(F_DMS, 0), '''51°30'49"N, 000°05'54"W''')
        d = LAX.destination(100, 66, radius=R_NM)
        self.test('destination', d.toStr(F_DM, prec=0), "34°37'N, 116°33'W")
        self.test('destination', d, '34.613643°N, 116.551171°W')

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            p = LatLon(53.2611, -0.7972)
            s = LatLon(53.3206, -1.7297)
            try:
                d = p.crossTrackDistanceTo(s, 96)
                self.test('crossTrackDistanceTo', d, '-305.67', '%.2f')  # -305.7
            except NotImplementedError as x:
                self.test('crossTrackDistanceTo', x, 'LatLon.crossTrackDistanceTo(end=bearing)')
            e = LatLon(53.1887, 0.1334)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', '%.2f')  # -307.5

        if hasattr(LatLon, 'greatCircle'):
            p = LatLon(53.3206, -1.7297)
            gc = p.greatCircle(96.0)
            self.test('greatCircle', gc, '(-0.79408, 0.12856, 0.59406)')

        if hasattr(LatLon, 'greatCircleTo'):
            p = LatLon(53.3206, -1.7297)
            q = LatLon(53.1887, 0.1334)
            gc = p.greatCircleTo(q)
            self.test('greatCircleTo', gc, '(-0.79408, 0.12859, 0.59406)')

    def testLatLonAttr(self, *modules):
        self.title('LatLon.attrs', __version__)
        attrs = {}
        for m in modules:
            ll = m.LatLon(0, 0)
            for a in dir(ll):
                if not a.startswith('__'):
                    a += _type(m.LatLon, a)
                    attrs[a] = attrs.get(a, ()) + (m.__name__,)
        for a, m in sorted(attrs.items()):
            m = ', '.join(sorted(m))
            self.test(a, m, m)  # passes always

        self.title('LatLon.mro', __version__)
        for m in modules:
            c = ', '.join(str(c)[8:-2] for c in m.LatLon.mro()[:-1])
            self.test(m.__name__, c, c)  # passes always

    def testModule(self, m, name=''):
        # check that __all__ names exist in module m
        self.title(m.__file__, m.__version__)
        for a in sorted(m.__all__):
            n = (name or m.__name__) + '.' + a + _type(m, a)
            t = getattr(m, a)
            o = getattr(t, '__module__', None)
            if o and o != m.__name__:
                n = '%s (%s)' % (n, o)
            self.test(n, hasattr(m, a), 'True')

    def testSpherical(self, LatLon, Nvector=None):
        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        i = p.intermediateTo(q, 0.25)
        self.test('intermediateTo', i, '51.372294°N, 000.707192°E' if Nvector
                                  else '51.372084°N, 000.707337°E')  # 51.3721°N, 000.7073°E  # Trig

        p = LatLon(51.8853, 0.2545)
        q = LatLon(49.0034, 2.5735)
        i = p.intersection(108.55, q, 32.44)
        self.test('intersection', i.toStr(F_D),  '50.907608°N, 004.508575°E' if Nvector
                                            else '50.907608°N, 004.508575°E')  # 50.9076°N, 004.5086°E  # Trig
        self.test('intersection', i.toStr(F_DMS), '50°54′27.39″N, 004°30′30.87″E')

        REO = LatLon(42.600, -117.866)
        BKE = LatLon(44.840, -117.806)
        i = REO.intersection(51, BKE, 137)
        self.test('intersection', i.toStr(F_D), '43.5719°N, 116.188757°W' if Nvector
                                           else '43.5719°N, 116.188757°W')  # 43.572°N, 116.189°W
        self.test('intersection', i.toStr(F_DMS), '43°34′18.84″N, 116°11′19.53″W')

        p = LatLon(0, 0)
        self.test('maxLat0',  p.maxLat( 0), '90.0')
        self.test('maxLat1',  p.maxLat( 1), '89.0')
        self.test('maxLat90', p.maxLat(90),  '0.0')

        if hasattr(LatLon, 'crossingParallels'):
            ps = p.crossingParallels(LatLon(60, 30), 30)
            t = ', '.join(map(lonDMS, ps))
            self.test('crossingParallels', t, '''009°35'38.65"E, 170°24'21.35"E''')

        p = LatLon(51.127, 1.338)
        q = LatLon(50.964, 1.853)
        b = p.rhumbBearingTo(q)
        self.test('rhumbBearingTo', b, '116.722', '%.3f')  # 116.7

        d = p.rhumbDistanceTo(q)
        self.test('rhumbDistanceTo', d, '40307.8', '%.1f')  # 40310 ?

        m = p.rhumbMidpointTo(q)
        self.test('rhumbMidpointo', m, '51.0455°N, 001.595727°E')  # 51.0455°N, 001.5957°E

    def testVectorial(self, LatLon, Nvector):
        if hasattr(LatLon, 'crossTrackDistanceTo'):
            p = LatLon(53.2611, -0.7972)
            s = LatLon(53.3206, -1.7297)
            d = p.crossTrackDistanceTo(s, 96.0)
            self.test('crossTrackDistanceTo', d, '-305.67', '%.2f')  # -305.7
            e = LatLon(53.1887, 0.1334)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', '%.2f')  # -307.5

        if hasattr(LatLon, 'enclosedby'):
            r = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            p = LatLon(45.1, 1.1)
            self.test('enclosedby', p.enclosedby(r), 'True')  # True
            r = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
            self.test('enclosedby', p.enclosedby(r), 'False')  # False

#       p = meanOf(r)
#       self.test('meanOf', p, '')

        v = Nvector(0.500, 0.500, 0.707)
        p = v.toLatLon()
        self.test('toLatLon', p, '44.995674°N, 045.0°E')  # 45.0°N, 45.0°E
        c = p.toNvector()
        self.test('toNvector', c, '(0.50004, 0.50004, 0.70705)')  # 0.500, 0.500, 0.707
        self.test('equals', c.equals(v), 'False')
        self.test('equals', c.equals(v, units=True), 'True')

        v = Nvector(52.205, 0.119, 0.0)
        c = v.copy()
        self.test('copy', c.equals(v), 'True')

    def testVincenty(self, LatLon, datum, VincentyError):
        d = datum
        n = ' (%s)' % (d.name,)

        Newport_RI = LatLon(41.49008, -71.312796, datum=d)
        Cleveland_OH = LatLon(41.499498, -81.695391, datum=d)
        m = Newport_RI.distanceTo(Cleveland_OH)
        self.test('distanceTo' + n, '%.5f' % m, '866455.43292')

        try:
            t = None
            m = Newport_RI.distanceTo(Newport_RI)
        except VincentyError as x:
            t = x  # Python 3+
        self.test('VincentyError' + n, t, 'LatLon(41.49008°N, 071.312796°W) coincident with LatLon(41.49008°N, 071.312796°W)')

        if hasattr(LatLon, 'toCartesian'):
            try:
                m = Newport_RI.distanceTo(Cleveland_OH.convertDatum(Datums.OSGB36))
                self.test('ValueError' + n, None, 'other Ellipsoid mistmatch: ...' + d.ellipsoid.name)
            except ValueError as x:
                self.test('ValueError' + n, x, 'other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.' + d.ellipsoid.name)
            except Exception as x:
                self.test('ValueError' + n, x, 'ValueError ...' + d.ellipsoid.name)

        p = LatLon(50.06632, -5.71475, datum=d)
        q = LatLon(58.64402, -3.07009, datum=d)
        m = p.distanceTo(q)
        self.test('distanceTo' + n, '%.4f' % m, '969954.1663')

        self.test('copy', p.copy().equals(p), 'True')

        t = p.distanceTo3(q)
        t = fStr(t, prec=6)
        self.test('distanceTo3' + n, t, '969954.166314, 9.141877, 11.29722')

        p = LatLon(37.95103, 144.42487, datum=d)
        q = LatLon(37.65280, 143.9265, datum=d)
        m = p.distanceTo(q)
        self.test('distanceTo' + n, '%.3f' % m, '54973.295')

        t = p.distanceTo3(q)
        t = fStr(t, prec=5)
        self.test('distanceTo3' + n, t, '54973.29527, 126.86992, 127.17539')

        p = LatLon(-37.95103, 144.42487, datum=d)
        p, f = p.destination2(54972.271, 306.86816)
        t = p.toStr(F_D) + ', ' + compassDMS(f, prec=4)
        self.test('destination2' + n, t, '37.652818°S, 143.926498°E, 307.1736°NW')


if __name__ == '__main__':

    import __init__ as init, datum, dms, utils, \
             ellipsoidalNvector, ellipsoidalVincenty, \
             sphericalNvector, sphericalTrigonometry, \
             nvector, vector3d   # PYCHOK expected

    t = Tests(__file__, __version__)
    # check that __all__ names exist in each module
    t.testModule(init, 'geodesy')
    for m in (datum, dms, utils,
              ellipsoidalNvector, ellipsoidalVincenty,
              sphericalNvector, sphericalTrigonometry,
              nvector, vector3d):
        t.testModule(m)
    t.testLatLonAttr(ellipsoidalNvector, ellipsoidalVincenty,
                     sphericalNvector, sphericalTrigonometry)
    t.results(nl=1)

    # testing tests.py version 16.09.06

    # testing __init__.pyc version 16.09.05
    # test 1 geodesy.Datum() class (datum): True
    # test 2 geodesy.Datums attribute (datum): True
    # test 3 geodesy.EPS float: True
    # test 4 geodesy.EPS1 float: True
    # test 5 geodesy.EPS2 float: True
    # test 6 geodesy.Ellipsoid() class (datum): True
    # test 7 geodesy.Ellipsoids attribute (datum): True
    # test 8 geodesy.F_D str: True
    # test 9 geodesy.F_DM str: True
    # test 10 geodesy.F_DMS str: True
    # test 11 geodesy.F_RAD str: True
    # test 12 geodesy.PI float: True
    # test 13 geodesy.PI2 float: True
    # test 14 geodesy.PI_2 float: True
    # test 15 geodesy.R_KM float: True
    # test 16 geodesy.R_M float: True
    # test 17 geodesy.R_NM float: True
    # test 18 geodesy.R_SM float: True
    # test 19 geodesy.S_DEG str: True
    # test 20 geodesy.S_MIN str: True
    # test 21 geodesy.S_SEC str: True
    # test 22 geodesy.S_SEP str: True
    # test 23 geodesy.Transform() class (datum): True
    # test 24 geodesy.Transforms attribute (datum): True
    # test 25 geodesy.VincentyError() class (ellipsoidalVincenty): True
    # test 26 geodesy.bearingDMS() function (dms): True
    # test 27 geodesy.cbrt() function (utils): True
    # test 28 geodesy.compassDMS() function (dms): True
    # test 29 geodesy.compassPoint() function (dms): True
    # test 30 geodesy.degrees attribute (math): True
    # test 31 geodesy.degrees180() function (utils): True
    # test 32 geodesy.degrees360() function (utils): True
    # test 33 geodesy.degrees90() function (utils): True
    # test 34 geodesy.ellipsoidalNvector module: True
    # test 35 geodesy.ellipsoidalVincenty module: True
    # test 36 geodesy.fStr() function (utils): True
    # test 37 geodesy.fdot() function (utils): True
    # test 38 geodesy.fsum attribute (math): True
    # test 39 geodesy.hypot3() function (utils): True
    # test 40 geodesy.isscalar() function (utils): True
    # test 41 geodesy.latDMS() function (dms): True
    # test 42 geodesy.len2() function (utils): True
    # test 43 geodesy.lonDMS() function (dms): True
    # test 44 geodesy.normDMS() function (dms): True
    # test 45 geodesy.parse3llh() function (dms): True
    # test 46 geodesy.parseDMS() function (dms): True
    # test 47 geodesy.precision() function (dms): True
    # test 48 geodesy.radians attribute (math): True
    # test 49 geodesy.radiansPI() function (utils): True
    # test 50 geodesy.radiansPI_2() function (utils): True
    # test 51 geodesy.sin_2() function (utils): True
    # test 52 geodesy.sphericalNvector module: True
    # test 53 geodesy.sphericalTrigonometry module: True
    # test 54 geodesy.tanPI_2_2() function (utils): True
    # test 55 geodesy.toDMS() function (dms): True
    # test 56 geodesy.wrapPI() function (utils): True
    # test 57 geodesy.wrapPI2() function (utils): True
    # test 58 geodesy.wrapPI_2() function (utils): True

    # testing datum.pyc version 16.09.05
    # test 59 datum.Datum() class: True
    # test 60 datum.Datums attribute: True
    # test 61 datum.Ellipsoid() class: True
    # test 62 datum.Ellipsoids attribute: True
    # test 63 datum.R_KM float: True
    # test 64 datum.R_M float: True
    # test 65 datum.R_NM float: True
    # test 66 datum.R_SM float: True
    # test 67 datum.Transform() class: True
    # test 68 datum.Transforms attribute: True

    # testing dms.pyc version 16.09.05
    # test 69 dms.F_D str: True
    # test 70 dms.F_DM str: True
    # test 71 dms.F_DMS str: True
    # test 72 dms.F_RAD str: True
    # test 73 dms.S_DEG str: True
    # test 74 dms.S_MIN str: True
    # test 75 dms.S_SEC str: True
    # test 76 dms.S_SEP str: True
    # test 77 dms.bearingDMS() function: True
    # test 78 dms.compassDMS() function: True
    # test 79 dms.compassPoint() function: True
    # test 80 dms.latDMS() function: True
    # test 81 dms.lonDMS() function: True
    # test 82 dms.normDMS() function: True
    # test 83 dms.parse3llh() function: True
    # test 84 dms.parseDMS() function: True
    # test 85 dms.precision() function: True
    # test 86 dms.toDMS() function: True

    # testing utils.pyc version 16.09.03
    # test 87 utils.EPS float: True
    # test 88 utils.EPS1 float: True
    # test 89 utils.EPS2 float: True
    # test 90 utils.PI float: True
    # test 91 utils.PI2 float: True
    # test 92 utils.PI_2 float: True
    # test 93 utils.cbrt() function: True
    # test 94 utils.degrees attribute (math): True
    # test 95 utils.degrees180() function: True
    # test 96 utils.degrees360() function: True
    # test 97 utils.degrees90() function: True
    # test 98 utils.fStr() function: True
    # test 99 utils.fdot() function: True
    # test 100 utils.fsum attribute (math): True
    # test 101 utils.hypot3() function: True
    # test 102 utils.isscalar() function: True
    # test 103 utils.len2() function: True
    # test 104 utils.radians attribute (math): True
    # test 105 utils.radiansPI() function: True
    # test 106 utils.radiansPI_2() function: True
    # test 107 utils.sin_2() function: True
    # test 108 utils.tanPI_2_2() function: True
    # test 109 utils.wrapPI() function: True
    # test 110 utils.wrapPI2() function: True
    # test 111 utils.wrapPI_2() function: True

    # testing ellipsoidalNvector.pyc version 16.09.06
    # test 112 ellipsoidalNvector.Cartesian() class: True
    # test 113 ellipsoidalNvector.LatLon() class: True
    # test 114 ellipsoidalNvector.Ned() class: True
    # test 115 ellipsoidalNvector.Nvector() class: True
    # test 116 ellipsoidalNvector.meanOf() function: True
    # test 117 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.pyc version 16.09.04
    # test 118 ellipsoidalVincenty.LatLon() class: True
    # test 119 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.pyc version 16.09.06
    # test 120 sphericalNvector.LatLon() class: True
    # test 121 sphericalNvector.areaOf() function: True
    # test 122 sphericalNvector.intersection() function: True
    # test 123 sphericalNvector.meanOf() function: True
    # test 124 sphericalNvector.triangulate() function: True
    # test 125 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.pyc version 16.09.06
    # test 126 sphericalTrigonometry.LatLon() class: True
    # test 127 sphericalTrigonometry.meanOf() function: True

    # testing nvector.pyc version 16.09.06
    # test 128 nvector.NorthPole attribute: True
    # test 129 nvector.Nvector() class: True
    # test 130 nvector.SouthPole attribute: True
    # test 131 nvector.sumOf() function: True

    # testing vector3d.pyc version 16.09.06
    # test 132 vector3d.Vector3d() class: True
    # test 133 vector3d.sumOf() function: True

    # testing LatLon.attrs version 16.09.06
    # test 134 Top() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 135 _Nv attribute: ellipsoidalNvector, sphericalNvector
    # test 136 _alter() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 137 _datum attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 138 _direct() method: ellipsoidalVincenty
    # test 139 _epsilon float: ellipsoidalVincenty
    # test 140 _gc3() method: sphericalNvector
    # test 141 _inverse() method: ellipsoidalVincenty
    # test 142 _iterations int: ellipsoidalVincenty
    # test 143 _r3 attribute: ellipsoidalNvector
    # test 144 _rhumb3() method: sphericalNvector, sphericalTrigonometry
    # test 145 _rotation3() method: ellipsoidalNvector
    # test 146 _v3d attribute: sphericalTrigonometry
    # test 147 bearingTo() method: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 148 convertDatum() method: ellipsoidalNvector, ellipsoidalVincenty
    # test 149 copy() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 150 crossTrackDistanceTo() method: sphericalNvector, sphericalTrigonometry
    # test 151 crossingParallels() method: sphericalTrigonometry
    # test 152 datum property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 153 deltaTo() method: ellipsoidalNvector
    # test 154 destination() method: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 155 destination2() method: ellipsoidalVincenty
    # test 156 destinationNed() method: ellipsoidalNvector
    # test 157 destinationPoint() method: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 158 distanceTo() method: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 159 distanceTo3() method: ellipsoidalVincenty
    # test 160 ellipsoid() method: ellipsoidalNvector, ellipsoidalVincenty
    # test 161 ellipsoids() method: ellipsoidalNvector, ellipsoidalVincenty
    # test 162 enclosedBy() method: sphericalNvector
    # test 163 epsilon property: ellipsoidalVincenty
    # test 164 equals() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 165 finalBearingOn() method: ellipsoidalVincenty
    # test 166 finalBearingTo() method: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 167 greatCircle() method: sphericalNvector, sphericalTrigonometry
    # test 168 greatCircleTo() method: sphericalNvector
    # test 169 height int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 170 initialBearingTo() method: ellipsoidalVincenty
    # test 171 intermediatePointTo() method: ellipsoidalNvector, sphericalNvector
    # test 172 intermediateTo() method: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 173 intersection() method: sphericalNvector, sphericalTrigonometry
    # test 174 isEnclosedBy() method: sphericalNvector, sphericalTrigonometry
    # test 175 isWithin() method: sphericalNvector, sphericalTrigonometry
    # test 176 isWithinExtent() method: sphericalNvector, sphericalTrigonometry
    # test 177 iterations property: ellipsoidalVincenty
    # test 178 lat attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 179 lon attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 180 maxLat() method: sphericalNvector, sphericalTrigonometry
    # test 181 midpointTo() method: sphericalNvector, sphericalTrigonometry
    # test 182 minLat() method: sphericalNvector, sphericalTrigonometry
    # test 183 nearestOn() method: sphericalNvector, sphericalTrigonometry
    # test 184 nearestPointOnSegment() method: sphericalNvector, sphericalTrigonometry
    # test 185 notImplemented() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 186 others() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 187 parse() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 188 rhumbBearingTo() method: sphericalNvector, sphericalTrigonometry
    # test 189 rhumbDistanceTo() method: sphericalNvector, sphericalTrigonometry
    # test 190 rhumbMidpointTo() method: sphericalNvector, sphericalTrigonometry
    # test 191 to3xyz() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 192 to4xyzh() method: ellipsoidalNvector, sphericalNvector
    # test 193 toCartesian() method: ellipsoidalNvector
    # test 194 toDatum() method: ellipsoidalNvector, ellipsoidalVincenty
    # test 195 toNvector() method: ellipsoidalNvector, sphericalNvector
    # test 196 toStr() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 197 toVector3d() method: sphericalTrigonometry
    # test 198 toradians() method: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 199 triangulate() method: sphericalNvector
    # test 200 trilaterate() method: sphericalNvector

    # testing LatLon.mro version 16.09.06
    # test 201 ellipsoidalNvector: ellipsoidalNvector.LatLon, nvector._LatLonNvectorBase, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 202 ellipsoidalVincenty: ellipsoidalVincenty.LatLon, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 203 sphericalNvector: sphericalNvector.LatLon, nvector._LatLonNvectorBase, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base
    # test 204 sphericalTrigonometry: sphericalTrigonometry.LatLon, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base

    # all tests.py tests passed (Python 2.7.10)

    # testing tests.py version 16.09.06

    # testing __init__.py version 16.09.05
    # test 1 geodesy.Datum() class (datum): True
    # test 2 geodesy.Datums attribute (datum): True
    # test 3 geodesy.EPS float: True
    # test 4 geodesy.EPS1 float: True
    # test 5 geodesy.EPS2 float: True
    # test 6 geodesy.Ellipsoid() class (datum): True
    # test 7 geodesy.Ellipsoids attribute (datum): True
    # test 8 geodesy.F_D str: True
    # test 9 geodesy.F_DM str: True
    # test 10 geodesy.F_DMS str: True
    # test 11 geodesy.F_RAD str: True
    # test 12 geodesy.PI float: True
    # test 13 geodesy.PI2 float: True
    # test 14 geodesy.PI_2 float: True
    # test 15 geodesy.R_KM float: True
    # test 16 geodesy.R_M float: True
    # test 17 geodesy.R_NM float: True
    # test 18 geodesy.R_SM float: True
    # test 19 geodesy.S_DEG str: True
    # test 20 geodesy.S_MIN str: True
    # test 21 geodesy.S_SEC str: True
    # test 22 geodesy.S_SEP str: True
    # test 23 geodesy.Transform() class (datum): True
    # test 24 geodesy.Transforms attribute (datum): True
    # test 25 geodesy.VincentyError() class (ellipsoidalVincenty): True
    # test 26 geodesy.bearingDMS() function (dms): True
    # test 27 geodesy.cbrt() function (utils): True
    # test 28 geodesy.compassDMS() function (dms): True
    # test 29 geodesy.compassPoint() function (dms): True
    # test 30 geodesy.degrees attribute (math): True
    # test 31 geodesy.degrees180() function (utils): True
    # test 32 geodesy.degrees360() function (utils): True
    # test 33 geodesy.degrees90() function (utils): True
    # test 34 geodesy.ellipsoidalNvector module: True
    # test 35 geodesy.ellipsoidalVincenty module: True
    # test 36 geodesy.fStr() function (utils): True
    # test 37 geodesy.fdot() function (utils): True
    # test 38 geodesy.fsum attribute (math): True
    # test 39 geodesy.hypot3() function (utils): True
    # test 40 geodesy.isscalar() function (utils): True
    # test 41 geodesy.latDMS() function (dms): True
    # test 42 geodesy.len2() function (utils): True
    # test 43 geodesy.lonDMS() function (dms): True
    # test 44 geodesy.normDMS() function (dms): True
    # test 45 geodesy.parse3llh() function (dms): True
    # test 46 geodesy.parseDMS() function (dms): True
    # test 47 geodesy.precision() function (dms): True
    # test 48 geodesy.radians attribute (math): True
    # test 49 geodesy.radiansPI() function (utils): True
    # test 50 geodesy.radiansPI_2() function (utils): True
    # test 51 geodesy.sin_2() function (utils): True
    # test 52 geodesy.sphericalNvector module: True
    # test 53 geodesy.sphericalTrigonometry module: True
    # test 54 geodesy.tanPI_2_2() function (utils): True
    # test 55 geodesy.toDMS() function (dms): True
    # test 56 geodesy.wrapPI() function (utils): True
    # test 57 geodesy.wrapPI2() function (utils): True
    # test 58 geodesy.wrapPI_2() function (utils): True

    # testing datum.py version 16.09.05
    # test 59 datum.Datum() class: True
    # test 60 datum.Datums attribute: True
    # test 61 datum.Ellipsoid() class: True
    # test 62 datum.Ellipsoids attribute: True
    # test 63 datum.R_KM float: True
    # test 64 datum.R_M float: True
    # test 65 datum.R_NM float: True
    # test 66 datum.R_SM float: True
    # test 67 datum.Transform() class: True
    # test 68 datum.Transforms attribute: True

    # testing dms.py version 16.09.05
    # test 69 dms.F_D str: True
    # test 70 dms.F_DM str: True
    # test 71 dms.F_DMS str: True
    # test 72 dms.F_RAD str: True
    # test 73 dms.S_DEG str: True
    # test 74 dms.S_MIN str: True
    # test 75 dms.S_SEC str: True
    # test 76 dms.S_SEP str: True
    # test 77 dms.bearingDMS() function: True
    # test 78 dms.compassDMS() function: True
    # test 79 dms.compassPoint() function: True
    # test 80 dms.latDMS() function: True
    # test 81 dms.lonDMS() function: True
    # test 82 dms.normDMS() function: True
    # test 83 dms.parse3llh() function: True
    # test 84 dms.parseDMS() function: True
    # test 85 dms.precision() function: True
    # test 86 dms.toDMS() function: True

    # testing utils.py version 16.09.03
    # test 87 utils.EPS float: True
    # test 88 utils.EPS1 float: True
    # test 89 utils.EPS2 float: True
    # test 90 utils.PI float: True
    # test 91 utils.PI2 float: True
    # test 92 utils.PI_2 float: True
    # test 93 utils.cbrt() function: True
    # test 94 utils.degrees attribute (math): True
    # test 95 utils.degrees180() function: True
    # test 96 utils.degrees360() function: True
    # test 97 utils.degrees90() function: True
    # test 98 utils.fStr() function: True
    # test 99 utils.fdot() function: True
    # test 100 utils.fsum attribute (math): True
    # test 101 utils.hypot3() function: True
    # test 102 utils.isscalar() function: True
    # test 103 utils.len2() function: True
    # test 104 utils.radians attribute (math): True
    # test 105 utils.radiansPI() function: True
    # test 106 utils.radiansPI_2() function: True
    # test 107 utils.sin_2() function: True
    # test 108 utils.tanPI_2_2() function: True
    # test 109 utils.wrapPI() function: True
    # test 110 utils.wrapPI2() function: True
    # test 111 utils.wrapPI_2() function: True

    # testing ellipsoidalNvector.py version 16.09.06
    # test 112 ellipsoidalNvector.Cartesian() class: True
    # test 113 ellipsoidalNvector.LatLon() class: True
    # test 114 ellipsoidalNvector.Ned() class: True
    # test 115 ellipsoidalNvector.Nvector() class: True
    # test 116 ellipsoidalNvector.meanOf() function: True
    # test 117 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.py version 16.09.04
    # test 118 ellipsoidalVincenty.LatLon() class: True
    # test 119 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.py version 16.09.06
    # test 120 sphericalNvector.LatLon() class: True
    # test 121 sphericalNvector.areaOf() function: True
    # test 122 sphericalNvector.intersection() function: True
    # test 123 sphericalNvector.meanOf() function: True
    # test 124 sphericalNvector.triangulate() function: True
    # test 125 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.py version 16.09.06
    # test 126 sphericalTrigonometry.LatLon() class: True
    # test 127 sphericalTrigonometry.meanOf() function: True

    # testing nvector.py version 16.09.06
    # test 128 nvector.NorthPole attribute: True
    # test 129 nvector.Nvector() class: True
    # test 130 nvector.SouthPole attribute: True
    # test 131 nvector.sumOf() function: True

    # testing vector3d.py version 16.09.06
    # test 132 vector3d.Vector3d() class: True
    # test 133 vector3d.sumOf() function: True

    # testing LatLon.attrs version 16.09.06
    # test 134 Top() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 135 _Nv attribute: ellipsoidalNvector, sphericalNvector
    # test 136 _alter() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 137 _datum attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 138 _direct() function: ellipsoidalVincenty
    # test 139 _epsilon float: ellipsoidalVincenty
    # test 140 _gc3() function: sphericalNvector
    # test 141 _inverse() function: ellipsoidalVincenty
    # test 142 _iterations int: ellipsoidalVincenty
    # test 143 _r3 attribute: ellipsoidalNvector
    # test 144 _rhumb3() function: sphericalNvector, sphericalTrigonometry
    # test 145 _rotation3() function: ellipsoidalNvector
    # test 146 _v3d attribute: sphericalTrigonometry
    # test 147 bearingTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 148 convertDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 149 copy() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 150 crossTrackDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 151 crossingParallels() function: sphericalTrigonometry
    # test 152 datum property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 153 deltaTo() function: ellipsoidalNvector
    # test 154 destination() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 155 destination2() function: ellipsoidalVincenty
    # test 156 destinationNed() function: ellipsoidalNvector
    # test 157 destinationPoint() function: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 158 distanceTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 159 distanceTo3() function: ellipsoidalVincenty
    # test 160 ellipsoid() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 161 ellipsoids() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 162 enclosedBy() function: sphericalNvector
    # test 163 epsilon property: ellipsoidalVincenty
    # test 164 equals() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 165 finalBearingOn() function: ellipsoidalVincenty
    # test 166 finalBearingTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 167 greatCircle() function: sphericalNvector, sphericalTrigonometry
    # test 168 greatCircleTo() function: sphericalNvector
    # test 169 height int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 170 initialBearingTo() function: ellipsoidalVincenty
    # test 171 intermediatePointTo() function: ellipsoidalNvector, sphericalNvector
    # test 172 intermediateTo() function: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 173 intersection() function: sphericalNvector, sphericalTrigonometry
    # test 174 isEnclosedBy() function: sphericalNvector, sphericalTrigonometry
    # test 175 isWithin() function: sphericalNvector, sphericalTrigonometry
    # test 176 isWithinExtent() function: sphericalNvector, sphericalTrigonometry
    # test 177 iterations property: ellipsoidalVincenty
    # test 178 lat attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 179 lon attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 180 maxLat() function: sphericalNvector, sphericalTrigonometry
    # test 181 midpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 182 minLat() function: sphericalNvector, sphericalTrigonometry
    # test 183 nearestOn() function: sphericalNvector, sphericalTrigonometry
    # test 184 nearestPointOnSegment() function: sphericalNvector, sphericalTrigonometry
    # test 185 notImplemented() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 186 others() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 187 parse() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 188 rhumbBearingTo() function: sphericalNvector, sphericalTrigonometry
    # test 189 rhumbDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 190 rhumbMidpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 191 to3xyz() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 192 to4xyzh() function: ellipsoidalNvector, sphericalNvector
    # test 193 toCartesian() function: ellipsoidalNvector
    # test 194 toDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 195 toNvector() function: ellipsoidalNvector, sphericalNvector
    # test 196 toStr() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 197 toVector3d() function: sphericalTrigonometry
    # test 198 toradians() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 199 triangulate() function: sphericalNvector
    # test 200 trilaterate() function: sphericalNvector

    # testing LatLon.mro version 16.09.06
    # test 201 ellipsoidalNvector: ellipsoidalNvector.LatLon, nvector._LatLonNvectorBase, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 202 ellipsoidalVincenty: ellipsoidalVincenty.LatLon, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 203 sphericalNvector: sphericalNvector.LatLon, nvector._LatLonNvectorBase, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base
    # test 204 sphericalTrigonometry: sphericalTrigonometry.LatLon, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base

    # all tests.py tests passed (Python 3.5.1)
