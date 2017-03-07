
# -*- coding: utf-8 -*-

# Tests for various geodesy modules.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# and <http://www.movable-type.co.uk/scripts/latlong.html>.

from os.path import basename, dirname
from platform import architecture
import sys
try:
    import geodesy as _  # PYCHOK expected
except ImportError:
    # extend sys.path to ../.. directory
    sys.path.insert(0, dirname(dirname(__file__)))
from inspect import isclass, isfunction, ismethod, ismodule

from geodesy import R_NM, F_D, F_DM, F_DMS, F_RAD, \
                    degrees, normDMS  # PYCHOK expected

__all__ = ('Tests',)
__version__ = '17.03.07'

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
    elif ismethod(t):
        t = '() method'
    elif isfunction(t):
        t = '() function'
    elif isinstance(t, _int):
        t = ' int'
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
    known  = 0
    total  = 0

    def __init__(self, file, version, mod=None):
        self._file = file
        if mod:
            self._name = mod.__name__
            self.title(self._name, mod.__version__)
        else:
            self._name = basename(file)
            self.title(file, version)

    def errors(self):
        return self.failed - self.known  # new failures

    def exit(self, errors=0):
        sys.exit(min(errors + self.errors(), 7))

    def printf(self, fmt, *args, **kwds):  # nl=0
        nl = '\n' * kwds.get('nl', 0)
        print(nl + self._prefix + (fmt % args))

    def results(self, nl=0):
        v = 'Python %s %s' % (sys.version.split()[0], architecture()[0])
        n = self.failed
        if n:
            p = '' if n == 1 else 's'
            k = self.known or ''
            if k:
                k = ', incl. %s KNOWN' % (k,)
            r = '(%.1f%%) FAILED%s' % (100.0 * n / self.total, k)
        else:
            n, p, r = 'all', 's', 'passed'
        self.printf('%s %s test%s %s (%s)', n, self._name, p, r, v, nl=nl)

    def test(self, name, value, expect, fmt='%s', known=False):
        self.total += 1  # tests
        f, v = '', fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect):
            self.failed += 1  # failures
            f = '  FAILED'
            if known:  # failed before
                self.known += 1
                f += ', KNOWN'
            f = '%s, expected %s' % (f, expect)
        self.printf('test %d %s: %s%s', self.total, name, v, f)

    def title(self, module, version):
        self.printf('testing %s version %s', basename(module), version, nl=1)

    def testLatLon(self, LatLon, Sph=True):
        # basic LatLon class tests
        p = LatLon(52.20472, 0.14056)
        self.test('lat/lonDMS', p, '52.20472°N, 000.14056°E')  # 52.20472°N, 000.14056°E
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 3),  '''52°12.283'N, 000°08.434'E''')
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 4),  '''52°12.2832'N, 000°08.4336'E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 0), '''52°12'17"N, 000°08'26"E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 1), '''52°12'17.0"N, 000°08'26.0"E''')
        self.test('lat/lonDMS F_RAD', p.toStr(F_RAD, 6), '0.911144N, 0.002453E')
        q = LatLon(*map(degrees, p.to2ab()))
        self.test('equals', q.equals(p), 'True')

        LAX = LatLon(33.+57./60, -(118.+24./60))
        JFK = LatLon(degrees(0.709186), -degrees(1.287762))

        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        self.test('equals', p.equals(q), 'False')

        if hasattr(LatLon, 'bearingTo'):
            b = p.bearingTo(q)
            self.test('bearingTo', b, '156.1666', '%.4f')  # 156.2
            b = p.finalBearingTo(q)
            self.test('finalBearingTo', b, '157.8904', '%.4f')
            b = LAX.bearingTo(JFK)
            self.test('bearingTo', b, '65.8921', '%.4f')  # 66

        c = p.copy()
        self.test('copy', p.equals(c), 'True')

        if hasattr(LatLon, 'distanceTo'):
            d = p.distanceTo(q)
            self.test('distanceTo', d, '404279.720589' if Sph else '404607.805988', '%.6f')  # 404300
            d = q.distanceTo(p)
            self.test('distanceTo', d, '404279.720589' if Sph else '404607.805988', '%.6f')  # 404300
            d = LAX.distanceTo(JFK, radius=R_NM) if Sph else LAX.distanceTo(JFK)
            self.test('distanceTo', d, '2145' if Sph else '3981601', '%.0f')  # PYCHOK false?

        if hasattr(LatLon, 'midpointTo'):
            m = p.midpointTo(q)
            self.test('midpointTo', m, '50.536327°N, 001.274614°E')  # PYCHOK false?  # 50.5363°N, 001.2746°E

        if hasattr(LatLon, 'destination'):
            p = LatLon(51.4778, -0.0015)
            d = p.destination(7794, 300.7)
            self.test('destination', d, '51.513546°N, 000.098345°W' if Sph else '51.513526°N, 000.098038°W')  # 51.5135°N, 0.0983°W ???
            self.test('destination', d.toStr(F_DMS, 0), '51°30′49″N, 000°05′54″W' if Sph else '51°30′49″N, 000°05′53″W')
            d = LAX.destination(100, 66, radius=R_NM) if Sph else LAX.destination(100, 66)
            self.test('destination', d.toStr(F_DM, prec=0), "34°37'N, 116°33'W" if Sph else "33°57'N, 118°24'W")
            self.test('destination', d, '34.613643°N, 116.551171°W' if Sph else '33.950367°N, 118.399012°W')  # PYCHOK false?

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            p = LatLon(53.2611, -0.7972)
            s = LatLon(53.3206, -1.7297)
            try:
                d = p.crossTrackDistanceTo(s, 96)
                self.test('crossTrackDistanceTo', d, '-305.67', '%.2f')  # -305.7
            except TypeError as x:
                self.test('crossTrackDistanceTo', x, 'type(end) mismatch: int vs sphericalTrigonometry.LatLon')  # PYCHOK false?
            e = LatLon(53.1887, 0.1334)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', '%.2f')  # PYCHOK false?  # -307.5

        if hasattr(LatLon, 'greatCircle'):
            p = LatLon(53.3206, -1.7297)
            gc = p.greatCircle(96.0)
            self.test('greatCircle', gc, '(-0.79408, 0.12856, 0.59406)')  # PYCHOK false?

        if hasattr(LatLon, 'greatCircleTo'):
            p = LatLon(53.3206, -1.7297)
            q = LatLon(53.1887, 0.1334)
            gc = p.greatCircleTo(q)
            self.test('greatCircleTo', gc, '(-0.79408, 0.12859, 0.59406)')  # PYCHOK false?

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
        m_ = (name or m.__name__).split('.')[-1] + '.'
        for a in sorted(m.__all__):
            n = m_ + a + _type(m, a)
            t = getattr(m, a)
            o = getattr(t, '__module__', None)
            if o and o != m.__name__:
                n = '%s (%s)' % (n, o)
            self.test(n, hasattr(m, a), 'True')

    def testVectorial(self, LatLon, Nvector, sumOf, isclockwise=None):
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

        class Nv(Nvector):
            pass
        v = Nvector(52.205, 0.119, 0.0)
        s = sumOf((v, c), Vector=Nv, h=0)
        self.test('sumOf', s, '(52.70504, 0.61904, 0.70705)')
        self.test('sumOf', s.__class__.__name__, 'Nv')

        c = v.copy()
        self.test('copy', c.equals(v), 'True')

        if isclockwise:
            f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            self.test('isclockwise', isclockwise(f), 'False')
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            self.test('isclockwise', isclockwise(t), 'True')
            try:
                self.test('isclockwise', isclockwise(t[:2]), ValueError)
            except ValueError as x:
                self.test('isclockwise', x, 'too few points: 2')

        if hasattr(LatLon, 'nearestOn'):
            s1 = LatLon(51.0, 1.0)
            s2 = LatLon(51.0, 2.0)
            s = LatLon(51.0, 1.9)
            p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E
            self.test('nearestOn', p.toStr(F_D, prec=4), '51.0004°N, 001.9°E')
            d = p.distanceTo(s)  # 42.71 m
            self.test('distanceTo', d, '42.712', fmt='%.3f')
            s = LatLon(51.0, 2.1)
            p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E
            self.test('nearestOn', p.toStr(F_D), '51.0°N, 002.0°E')

            # courtesy AkimboEG on GitHub
            s1 = LatLon(0, 0)
            s2 = LatLon(0, 1)
            s = LatLon(1, 0)
            p = s.nearestOn(s1, s2)  # 0.0°N, 0.0°E
            self.test('nearestOn', p, '00.0°N, 000.0°E')

            p = LatLon(10, -140).nearestOn(LatLon(0, 20), LatLon(0, 40))
            self.test('nearestOn', p, '00.0°N, 020.0°E')

        if hasattr(LatLon, 'triangulate'):
            # courtesy of pvezid  Feb 10, 2017
            p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
            self.test('BasseC', p, '47.3038°N, 002.5721°W')
            s = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
            self.test('BasseH', s, '47.311067°N, 002.528617°W')
            t = p.triangulate(7, s, 295)
            self.test('triangulate', t, '47.323667°N, 002.568501°W')


if __name__ == '__main__':

    from geodesy import datum, dms, lcc, mgrs, osgr, \
                        ellipsoidalNvector, ellipsoidalVincenty, \
                        sphericalNvector, sphericalTrigonometry, \
                        nvector, vector3d, utm, utils  # PYCHOK expected
    import geodesy  # PYCHOK expected

    t = Tests(__file__, __version__)
    # check that __all__ names exist in each module
    t.testModule(geodesy, 'geodesy')
    for m in (datum, dms, lcc, mgrs, osgr,
              ellipsoidalNvector, ellipsoidalVincenty,
              sphericalNvector, sphericalTrigonometry,
              nvector, vector3d, utm, utils):
        t.testModule(m)
    t.testLatLonAttr(ellipsoidalNvector, ellipsoidalVincenty,
                     sphericalNvector, sphericalTrigonometry)
    t.results(nl=1)
    t.exit()

    # Typical test results (on MacOS 10.12.3)

    # testing tests.py version 17.03.07

    # testing __init__.py version 17.03.07
    # test 1 geodesy.Conic() class (geodesy.lcc): True
    # test 2 geodesy.Conics attribute (geodesy.datum): True
    # test 3 geodesy.Datum() class (geodesy.datum): True
    # test 4 geodesy.Datums attribute (geodesy.datum): True
    # test 5 geodesy.EPS float: True
    # test 6 geodesy.EPS1 float: True
    # test 7 geodesy.EPS2 float: True
    # test 8 geodesy.Ellipsoid() class (geodesy.datum): True
    # test 9 geodesy.Ellipsoids attribute (geodesy.datum): True
    # test 10 geodesy.F_D str: True
    # test 11 geodesy.F_DM str: True
    # test 12 geodesy.F_DMS str: True
    # test 13 geodesy.F_RAD str: True
    # test 14 geodesy.Lcc() class (geodesy.lcc): True
    # test 15 geodesy.Mgrs() class (geodesy.mgrs): True
    # test 16 geodesy.Osgr() class (geodesy.osgr): True
    # test 17 geodesy.PI float: True
    # test 18 geodesy.PI2 float: True
    # test 19 geodesy.PI_2 float: True
    # test 20 geodesy.R_KM float: True
    # test 21 geodesy.R_M float: True
    # test 22 geodesy.R_NM float: True
    # test 23 geodesy.R_SM float: True
    # test 24 geodesy.S_DEG str: True
    # test 25 geodesy.S_MIN str: True
    # test 26 geodesy.S_SEC str: True
    # test 27 geodesy.S_SEP str: True
    # test 28 geodesy.Transform() class (geodesy.datum): True
    # test 29 geodesy.Transforms attribute (geodesy.datum): True
    # test 30 geodesy.Utm() class (geodesy.utm): True
    # test 31 geodesy.VincentyError() class (geodesy.ellipsoidalVincenty): True
    # test 32 geodesy.bearingDMS() function (geodesy.dms): True
    # test 33 geodesy.cbrt() function (geodesy.utils): True
    # test 34 geodesy.compassDMS() function (geodesy.dms): True
    # test 35 geodesy.compassPoint() function (geodesy.dms): True
    # test 36 geodesy.degrees attribute (math): True
    # test 37 geodesy.degrees180() function (geodesy.utils): True
    # test 38 geodesy.degrees360() function (geodesy.utils): True
    # test 39 geodesy.degrees90() function (geodesy.utils): True
    # test 40 geodesy.ellipsoidalNvector module: True
    # test 41 geodesy.ellipsoidalVincenty module: True
    # test 42 geodesy.fStr() function (geodesy.utils): True
    # test 43 geodesy.false2f() function (geodesy.utils): True
    # test 44 geodesy.fdot() function (geodesy.utils): True
    # test 45 geodesy.fdot3() function (geodesy.utils): True
    # test 46 geodesy.fsum attribute (math): True
    # test 47 geodesy.halfs() function (geodesy.utils): True
    # test 48 geodesy.hsin() function (geodesy.utils): True
    # test 49 geodesy.hypot1() function (geodesy.utils): True
    # test 50 geodesy.hypot3() function (geodesy.utils): True
    # test 51 geodesy.isint() function (geodesy.utils): True
    # test 52 geodesy.isscalar() function (geodesy.utils): True
    # test 53 geodesy.latDMS() function (geodesy.dms): True
    # test 54 geodesy.len2() function (geodesy.utils): True
    # test 55 geodesy.lonDMS() function (geodesy.dms): True
    # test 56 geodesy.map2() function (geodesy.utils): True
    # test 57 geodesy.normDMS() function (geodesy.dms): True
    # test 58 geodesy.nvector module: True
    # test 59 geodesy.parse3llh() function (geodesy.dms): True
    # test 60 geodesy.parseDMS() function (geodesy.dms): True
    # test 61 geodesy.parseMGRS() function (geodesy.mgrs): True
    # test 62 geodesy.parseOSGR() function (geodesy.osgr): True
    # test 63 geodesy.parseUTM() function (geodesy.utm): True
    # test 64 geodesy.precision() function (geodesy.dms): True
    # test 65 geodesy.radians attribute (math): True
    # test 66 geodesy.radiansPI() function (geodesy.utils): True
    # test 67 geodesy.radiansPI2() function (geodesy.utils): True
    # test 68 geodesy.radiansPI_2() function (geodesy.utils): True
    # test 69 geodesy.sphericalNvector module: True
    # test 70 geodesy.sphericalTrigonometry module: True
    # test 71 geodesy.tanPI_2_2() function (geodesy.utils): True
    # test 72 geodesy.toDMS() function (geodesy.dms): True
    # test 73 geodesy.toLcc() function (geodesy.lcc): True
    # test 74 geodesy.toMgrs() function (geodesy.mgrs): True
    # test 75 geodesy.toOsgr() function (geodesy.osgr): True
    # test 76 geodesy.toUtm() function (geodesy.utm): True
    # test 77 geodesy.vector3d module: True
    # test 78 geodesy.wrap180() function (geodesy.utils): True
    # test 79 geodesy.wrap90() function (geodesy.utils): True
    # test 80 geodesy.wrapPI() function (geodesy.utils): True
    # test 81 geodesy.wrapPI2() function (geodesy.utils): True
    # test 82 geodesy.wrapPI_2() function (geodesy.utils): True

    # testing datum.py version 17.02.27
    # test 83 datum.Datum() class: True
    # test 84 datum.Datums attribute: True
    # test 85 datum.Ellipsoid() class: True
    # test 86 datum.Ellipsoids attribute: True
    # test 87 datum.R_KM float: True
    # test 88 datum.R_M float: True
    # test 89 datum.R_NM float: True
    # test 90 datum.R_SM float: True
    # test 91 datum.Transform() class: True
    # test 92 datum.Transforms attribute: True

    # testing dms.py version 17.02.15
    # test 93 dms.F_D str: True
    # test 94 dms.F_DM str: True
    # test 95 dms.F_DMS str: True
    # test 96 dms.F_RAD str: True
    # test 97 dms.S_DEG str: True
    # test 98 dms.S_MIN str: True
    # test 99 dms.S_SEC str: True
    # test 100 dms.S_SEP str: True
    # test 101 dms.bearingDMS() function: True
    # test 102 dms.compassDMS() function: True
    # test 103 dms.compassPoint() function: True
    # test 104 dms.latDMS() function: True
    # test 105 dms.lonDMS() function: True
    # test 106 dms.normDMS() function: True
    # test 107 dms.parse3llh() function: True
    # test 108 dms.parseDMS() function: True
    # test 109 dms.precision() function: True
    # test 110 dms.toDMS() function: True

    # testing lcc.py version 17.02.14
    # test 111 lcc.Conic() class: True
    # test 112 lcc.Conics attribute (geodesy.datum): True
    # test 113 lcc.Lcc() class: True
    # test 114 lcc.toLcc() function: True

    # testing mgrs.py version 17.03.07
    # test 115 mgrs.Mgrs() class: True
    # test 116 mgrs.parseMGRS() function: True
    # test 117 mgrs.toMgrs() function: True

    # testing osgr.py version 17.03.07
    # test 118 osgr.Osgr() class: True
    # test 119 osgr.parseOSGR() function: True
    # test 120 osgr.toOsgr() function: True

    # testing ellipsoidalNvector.py version 17.02.15
    # test 121 ellipsoidalNvector.Cartesian() class: True
    # test 122 ellipsoidalNvector.LatLon() class: True
    # test 123 ellipsoidalNvector.Ned() class: True
    # test 124 ellipsoidalNvector.Nvector() class: True
    # test 125 ellipsoidalNvector.meanOf() function: True
    # test 126 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.py version 17.02.14
    # test 127 ellipsoidalVincenty.Cartesian() class: True
    # test 128 ellipsoidalVincenty.LatLon() class: True
    # test 129 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.py version 17.02.27
    # test 130 sphericalNvector.LatLon() class: True
    # test 131 sphericalNvector.Nvector() class: True
    # test 132 sphericalNvector.areaOf() function: True
    # test 133 sphericalNvector.intersection() function: True
    # test 134 sphericalNvector.isclockwise() function: True
    # test 135 sphericalNvector.meanOf() function: True
    # test 136 sphericalNvector.triangulate() function: True
    # test 137 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.py version 17.03.07
    # test 138 sphericalTrigonometry.LatLon() class: True
    # test 139 sphericalTrigonometry.intersection() function: True
    # test 140 sphericalTrigonometry.meanOf() function: True

    # testing nvector.py version 17.02.14
    # test 141 nvector.NorthPole attribute: True
    # test 142 nvector.Nvector() class: True
    # test 143 nvector.SouthPole attribute: True
    # test 144 nvector.sumOf() function: True

    # testing vector3d.py version 17.02.14
    # test 145 vector3d.Vector3d() class: True
    # test 146 vector3d.sumOf() function: True

    # testing utm.py version 17.03.07
    # test 147 utm.Utm() class: True
    # test 148 utm.parseUTM() function: True
    # test 149 utm.toUtm() function: True

    # testing utils.py version 17.03.07
    # test 150 utils.EPS float: True
    # test 151 utils.EPS1 float: True
    # test 152 utils.EPS2 float: True
    # test 153 utils.PI float: True
    # test 154 utils.PI2 float: True
    # test 155 utils.PI_2 float: True
    # test 156 utils.cbrt() function: True
    # test 157 utils.degrees attribute (math): True
    # test 158 utils.degrees180() function: True
    # test 159 utils.degrees360() function: True
    # test 160 utils.degrees90() function: True
    # test 161 utils.fStr() function: True
    # test 162 utils.false2f() function: True
    # test 163 utils.fdot() function: True
    # test 164 utils.fdot3() function: True
    # test 165 utils.fsum attribute (math): True
    # test 166 utils.halfs() function: True
    # test 167 utils.hsin() function: True
    # test 168 utils.hypot1() function: True
    # test 169 utils.hypot3() function: True
    # test 170 utils.isint() function: True
    # test 171 utils.isscalar() function: True
    # test 172 utils.len2() function: True
    # test 173 utils.map2() function: True
    # test 174 utils.radians attribute (math): True
    # test 175 utils.radiansPI() function: True
    # test 176 utils.radiansPI2() function: True
    # test 177 utils.radiansPI_2() function: True
    # test 178 utils.tanPI_2_2() function: True
    # test 179 utils.wrap180() function: True
    # test 180 utils.wrap90() function: True
    # test 181 utils.wrapPI() function: True
    # test 182 utils.wrapPI2() function: True
    # test 183 utils.wrapPI_2() function: True

    # testing LatLon.attrs version 17.03.07
    # test 184 _Nv attribute: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 185 _ab attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 186 _alter() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 187 _datum attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 188 _direct() method: geodesy.ellipsoidalVincenty
    # test 189 _epsilon float: geodesy.ellipsoidalVincenty
    # test 190 _gc3() method: geodesy.sphericalNvector
    # test 191 _height int: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 192 _inverse() method: geodesy.ellipsoidalVincenty
    # test 193 _iterations int: geodesy.ellipsoidalVincenty
    # test 194 _lat int: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 195 _lon int: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 196 _osgr attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 197 _r3 attribute: geodesy.ellipsoidalNvector
    # test 198 _rhumb3() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 199 _rotation3() method: geodesy.ellipsoidalNvector
    # test 200 _update() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 201 _utm attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 202 _v3d attribute: geodesy.sphericalTrigonometry
    # test 203 bearingTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 204 classname() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 205 convertDatum() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 206 copy() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 207 crossTrackDistanceTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 208 crossingParallels() method: geodesy.sphericalTrigonometry
    # test 209 datum property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 210 deltaTo() method: geodesy.ellipsoidalNvector
    # test 211 destination() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 212 destination2() method: geodesy.ellipsoidalVincenty
    # test 213 destinationNed() method: geodesy.ellipsoidalNvector
    # test 214 distanceTo() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 215 distanceTo3() method: geodesy.ellipsoidalVincenty
    # test 216 ellipsoid() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 217 ellipsoids() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 218 epsilon property: geodesy.ellipsoidalVincenty
    # test 219 equals() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 220 finalBearingOn() method: geodesy.ellipsoidalVincenty
    # test 221 finalBearingTo() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 222 greatCircle() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 223 greatCircleTo() method: geodesy.sphericalNvector
    # test 224 height property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 225 initialBearingTo() method: geodesy.ellipsoidalVincenty
    # test 226 intermediateChordTo() method: geodesy.sphericalNvector
    # test 227 intermediateTo() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 228 intersection() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 229 isEnclosedBy() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 230 isWithin() method: geodesy.sphericalNvector
    # test 231 isclockwise() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 232 iterations property: geodesy.ellipsoidalVincenty
    # test 233 lat property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 234 lon property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 235 maxLat() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 236 midpointTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 237 minLat() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 238 nearestOn() method: geodesy.sphericalNvector
    # test 239 others() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 240 parse() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 241 points() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 242 rhumbBearingTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 243 rhumbDistanceTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 244 rhumbMidpointTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 245 to2ab() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 246 to3xyz() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 247 to4xyzh() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 248 toCartesian() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 249 toDatum() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 250 toNvector() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 251 toOsgr() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 252 toStr() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 253 toStr2() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 254 toUtm() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 255 toVector3d() method: geodesy.sphericalTrigonometry
    # test 256 topsub() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 257 triangulate() method: geodesy.sphericalNvector
    # test 258 trilaterate() method: geodesy.sphericalNvector

    # testing LatLon.mro version 17.03.07
    # test 259 geodesy.ellipsoidalNvector: geodesy.ellipsoidalNvector.LatLon, geodesy.nvector.LatLonNvectorBase, geodesy.ellipsoidalBase.LatLonEllipsoidalBase, geodesy.bases.LatLonHeightBase, geodesy.bases.Base
    # test 260 geodesy.ellipsoidalVincenty: geodesy.ellipsoidalVincenty.LatLon, geodesy.ellipsoidalBase.LatLonEllipsoidalBase, geodesy.bases.LatLonHeightBase, geodesy.bases.Base
    # test 261 geodesy.sphericalNvector: geodesy.sphericalNvector.LatLon, geodesy.nvector.LatLonNvectorBase, geodesy.sphericalBase.LatLonSphericalBase, geodesy.bases.LatLonHeightBase, geodesy.bases.Base
    # test 262 geodesy.sphericalTrigonometry: geodesy.sphericalTrigonometry.LatLon, geodesy.sphericalBase.LatLonSphericalBase, geodesy.bases.LatLonHeightBase, geodesy.bases.Base

    # all tests.py tests passed (Python 2.7.13 64bit)

    # testing tests.py version 17.03.07

    # testing __init__.py version 17.03.07
    # test 1 geodesy.Conic() class (lcc): True
    # test 2 geodesy.Conics attribute (datum): True
    # test 3 geodesy.Datum() class (datum): True
    # test 4 geodesy.Datums attribute (datum): True
    # test 5 geodesy.EPS float: True
    # test 6 geodesy.EPS1 float: True
    # test 7 geodesy.EPS2 float: True
    # test 8 geodesy.Ellipsoid() class (datum): True
    # test 9 geodesy.Ellipsoids attribute (datum): True
    # test 10 geodesy.F_D str: True
    # test 11 geodesy.F_DM str: True
    # test 12 geodesy.F_DMS str: True
    # test 13 geodesy.F_RAD str: True
    # test 14 geodesy.Lcc() class (lcc): True
    # test 15 geodesy.Mgrs() class (mgrs): True
    # test 16 geodesy.Osgr() class (osgr): True
    # test 17 geodesy.PI float: True
    # test 18 geodesy.PI2 float: True
    # test 19 geodesy.PI_2 float: True
    # test 20 geodesy.R_KM float: True
    # test 21 geodesy.R_M float: True
    # test 22 geodesy.R_NM float: True
    # test 23 geodesy.R_SM float: True
    # test 24 geodesy.S_DEG str: True
    # test 25 geodesy.S_MIN str: True
    # test 26 geodesy.S_SEC str: True
    # test 27 geodesy.S_SEP str: True
    # test 28 geodesy.Transform() class (datum): True
    # test 29 geodesy.Transforms attribute (datum): True
    # test 30 geodesy.Utm() class (utm): True
    # test 31 geodesy.VincentyError() class (ellipsoidalVincenty): True
    # test 32 geodesy.bearingDMS() function (dms): True
    # test 33 geodesy.cbrt() function (utils): True
    # test 34 geodesy.compassDMS() function (dms): True
    # test 35 geodesy.compassPoint() function (dms): True
    # test 36 geodesy.degrees attribute (math): True
    # test 37 geodesy.degrees180() function (utils): True
    # test 38 geodesy.degrees360() function (utils): True
    # test 39 geodesy.degrees90() function (utils): True
    # test 40 geodesy.ellipsoidalNvector module: True
    # test 41 geodesy.ellipsoidalVincenty module: True
    # test 42 geodesy.fStr() function (utils): True
    # test 43 geodesy.false2f() function (utils): True
    # test 44 geodesy.fdot() function (utils): True
    # test 45 geodesy.fdot3() function (utils): True
    # test 46 geodesy.fsum attribute (math): True
    # test 47 geodesy.halfs() function (utils): True
    # test 48 geodesy.hsin() function (utils): True
    # test 49 geodesy.hypot1() function (utils): True
    # test 50 geodesy.hypot3() function (utils): True
    # test 51 geodesy.isint() function (utils): True
    # test 52 geodesy.isscalar() function (utils): True
    # test 53 geodesy.latDMS() function (dms): True
    # test 54 geodesy.len2() function (utils): True
    # test 55 geodesy.lonDMS() function (dms): True
    # test 56 geodesy.map2() function (utils): True
    # test 57 geodesy.normDMS() function (dms): True
    # test 58 geodesy.nvector module: True
    # test 59 geodesy.parse3llh() function (dms): True
    # test 60 geodesy.parseDMS() function (dms): True
    # test 61 geodesy.parseMGRS() function (mgrs): True
    # test 62 geodesy.parseOSGR() function (osgr): True
    # test 63 geodesy.parseUTM() function (utm): True
    # test 64 geodesy.precision() function (dms): True
    # test 65 geodesy.radians attribute (math): True
    # test 66 geodesy.radiansPI() function (utils): True
    # test 67 geodesy.radiansPI2() function (utils): True
    # test 68 geodesy.radiansPI_2() function (utils): True
    # test 69 geodesy.sphericalNvector module: True
    # test 70 geodesy.sphericalTrigonometry module: True
    # test 71 geodesy.tanPI_2_2() function (utils): True
    # test 72 geodesy.toDMS() function (dms): True
    # test 73 geodesy.toLcc() function (lcc): True
    # test 74 geodesy.toMgrs() function (mgrs): True
    # test 75 geodesy.toOsgr() function (osgr): True
    # test 76 geodesy.toUtm() function (utm): True
    # test 77 geodesy.vector3d module: True
    # test 78 geodesy.wrap180() function (utils): True
    # test 79 geodesy.wrap90() function (utils): True
    # test 80 geodesy.wrapPI() function (utils): True
    # test 81 geodesy.wrapPI2() function (utils): True
    # test 82 geodesy.wrapPI_2() function (utils): True

    # testing datum.py version 17.02.27
    # test 83 datum.Datum() class: True
    # test 84 datum.Datums attribute: True
    # test 85 datum.Ellipsoid() class: True
    # test 86 datum.Ellipsoids attribute: True
    # test 87 datum.R_KM float: True
    # test 88 datum.R_M float: True
    # test 89 datum.R_NM float: True
    # test 90 datum.R_SM float: True
    # test 91 datum.Transform() class: True
    # test 92 datum.Transforms attribute: True

    # testing dms.py version 17.02.15
    # test 93 dms.F_D str: True
    # test 94 dms.F_DM str: True
    # test 95 dms.F_DMS str: True
    # test 96 dms.F_RAD str: True
    # test 97 dms.S_DEG str: True
    # test 98 dms.S_MIN str: True
    # test 99 dms.S_SEC str: True
    # test 100 dms.S_SEP str: True
    # test 101 dms.bearingDMS() function: True
    # test 102 dms.compassDMS() function: True
    # test 103 dms.compassPoint() function: True
    # test 104 dms.latDMS() function: True
    # test 105 dms.lonDMS() function: True
    # test 106 dms.normDMS() function: True
    # test 107 dms.parse3llh() function: True
    # test 108 dms.parseDMS() function: True
    # test 109 dms.precision() function: True
    # test 110 dms.toDMS() function: True

    # testing lcc.py version 17.02.14
    # test 111 lcc.Conic() class: True
    # test 112 lcc.Conics attribute (datum): True
    # test 113 lcc.Lcc() class: True
    # test 114 lcc.toLcc() function: True

    # testing mgrs.py version 17.03.07
    # test 115 mgrs.Mgrs() class: True
    # test 116 mgrs.parseMGRS() function: True
    # test 117 mgrs.toMgrs() function: True

    # testing osgr.py version 17.03.07
    # test 118 osgr.Osgr() class: True
    # test 119 osgr.parseOSGR() function: True
    # test 120 osgr.toOsgr() function: True

    # testing ellipsoidalNvector.py version 17.02.15
    # test 121 ellipsoidalNvector.Cartesian() class: True
    # test 122 ellipsoidalNvector.LatLon() class: True
    # test 123 ellipsoidalNvector.Ned() class: True
    # test 124 ellipsoidalNvector.Nvector() class: True
    # test 125 ellipsoidalNvector.meanOf() function: True
    # test 126 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.py version 17.02.14
    # test 127 ellipsoidalVincenty.Cartesian() class: True
    # test 128 ellipsoidalVincenty.LatLon() class: True
    # test 129 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.py version 17.02.27
    # test 130 sphericalNvector.LatLon() class: True
    # test 131 sphericalNvector.Nvector() class: True
    # test 132 sphericalNvector.areaOf() function: True
    # test 133 sphericalNvector.intersection() function: True
    # test 134 sphericalNvector.isclockwise() function: True
    # test 135 sphericalNvector.meanOf() function: True
    # test 136 sphericalNvector.triangulate() function: True
    # test 137 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.py version 17.03.07
    # test 138 sphericalTrigonometry.LatLon() class: True
    # test 139 sphericalTrigonometry.intersection() function: True
    # test 140 sphericalTrigonometry.meanOf() function: True

    # testing nvector.py version 17.02.14
    # test 141 nvector.NorthPole attribute: True
    # test 142 nvector.Nvector() class: True
    # test 143 nvector.SouthPole attribute: True
    # test 144 nvector.sumOf() function: True

    # testing vector3d.py version 17.02.14
    # test 145 vector3d.Vector3d() class: True
    # test 146 vector3d.sumOf() function: True

    # testing utm.py version 17.03.07
    # test 147 utm.Utm() class: True
    # test 148 utm.parseUTM() function: True
    # test 149 utm.toUtm() function: True

    # testing utils.py version 17.03.07
    # test 150 utils.EPS float: True
    # test 151 utils.EPS1 float: True
    # test 152 utils.EPS2 float: True
    # test 153 utils.PI float: True
    # test 154 utils.PI2 float: True
    # test 155 utils.PI_2 float: True
    # test 156 utils.cbrt() function: True
    # test 157 utils.degrees attribute (math): True
    # test 158 utils.degrees180() function: True
    # test 159 utils.degrees360() function: True
    # test 160 utils.degrees90() function: True
    # test 161 utils.fStr() function: True
    # test 162 utils.false2f() function: True
    # test 163 utils.fdot() function: True
    # test 164 utils.fdot3() function: True
    # test 165 utils.fsum attribute (math): True
    # test 166 utils.halfs() function: True
    # test 167 utils.hsin() function: True
    # test 168 utils.hypot1() function: True
    # test 169 utils.hypot3() function: True
    # test 170 utils.isint() function: True
    # test 171 utils.isscalar() function: True
    # test 172 utils.len2() function: True
    # test 173 utils.map2() function: True
    # test 174 utils.radians attribute (math): True
    # test 175 utils.radiansPI() function: True
    # test 176 utils.radiansPI2() function: True
    # test 177 utils.radiansPI_2() function: True
    # test 178 utils.tanPI_2_2() function: True
    # test 179 utils.wrap180() function: True
    # test 180 utils.wrap90() function: True
    # test 181 utils.wrapPI() function: True
    # test 182 utils.wrapPI2() function: True
    # test 183 utils.wrapPI_2() function: True

    # testing LatLon.attrs version 17.03.07
    # test 184 _Nv attribute: ellipsoidalNvector, sphericalNvector
    # test 185 _ab attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 186 _alter() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 187 _datum attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 188 _direct() function: ellipsoidalVincenty
    # test 189 _epsilon float: ellipsoidalVincenty
    # test 190 _gc3() function: sphericalNvector
    # test 191 _height int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 192 _inverse() function: ellipsoidalVincenty
    # test 193 _iterations int: ellipsoidalVincenty
    # test 194 _lat int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 195 _lon int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 196 _osgr attribute: ellipsoidalNvector, ellipsoidalVincenty
    # test 197 _r3 attribute: ellipsoidalNvector
    # test 198 _rhumb3() function: sphericalNvector, sphericalTrigonometry
    # test 199 _rotation3() function: ellipsoidalNvector
    # test 200 _update() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 201 _utm attribute: ellipsoidalNvector, ellipsoidalVincenty
    # test 202 _v3d attribute: sphericalTrigonometry
    # test 203 bearingTo() function: sphericalNvector, sphericalTrigonometry
    # test 204 classname() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 205 convertDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 206 copy() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 207 crossTrackDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 208 crossingParallels() function: sphericalTrigonometry
    # test 209 datum property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 210 deltaTo() function: ellipsoidalNvector
    # test 211 destination() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 212 destination2() function: ellipsoidalVincenty
    # test 213 destinationNed() function: ellipsoidalNvector
    # test 214 distanceTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 215 distanceTo3() function: ellipsoidalVincenty
    # test 216 ellipsoid() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 217 ellipsoids() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 218 epsilon property: ellipsoidalVincenty
    # test 219 equals() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 220 finalBearingOn() function: ellipsoidalVincenty
    # test 221 finalBearingTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 222 greatCircle() function: sphericalNvector, sphericalTrigonometry
    # test 223 greatCircleTo() function: sphericalNvector
    # test 224 height property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 225 initialBearingTo() function: ellipsoidalVincenty
    # test 226 intermediateChordTo() function: sphericalNvector
    # test 227 intermediateTo() function: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 228 intersection() function: sphericalNvector, sphericalTrigonometry
    # test 229 isEnclosedBy() function: sphericalNvector, sphericalTrigonometry
    # test 230 isWithin() function: sphericalNvector
    # test 231 isclockwise() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 232 iterations property: ellipsoidalVincenty
    # test 233 lat property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 234 lon property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 235 maxLat() function: sphericalNvector, sphericalTrigonometry
    # test 236 midpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 237 minLat() function: sphericalNvector, sphericalTrigonometry
    # test 238 nearestOn() function: sphericalNvector
    # test 239 others() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 240 parse() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 241 points() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 242 rhumbBearingTo() function: sphericalNvector, sphericalTrigonometry
    # test 243 rhumbDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 244 rhumbMidpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 245 to2ab() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 246 to3xyz() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 247 to4xyzh() function: ellipsoidalNvector, sphericalNvector
    # test 248 toCartesian() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 249 toDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 250 toNvector() function: ellipsoidalNvector, sphericalNvector
    # test 251 toOsgr() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 252 toStr() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 253 toStr2() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 254 toUtm() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 255 toVector3d() function: sphericalTrigonometry
    # test 256 topsub() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 257 triangulate() function: sphericalNvector
    # test 258 trilaterate() function: sphericalNvector

    # testing LatLon.mro version 17.03.07
    # test 259 ellipsoidalNvector: ellipsoidalNvector.LatLon, nvector.LatLonNvectorBase, ellipsoidalBase.LatLonEllipsoidalBase, bases.LatLonHeightBase, bases.Base
    # test 260 ellipsoidalVincenty: ellipsoidalVincenty.LatLon, ellipsoidalBase.LatLonEllipsoidalBase, bases.LatLonHeightBase, bases.Base
    # test 261 sphericalNvector: sphericalNvector.LatLon, nvector.LatLonNvectorBase, sphericalBase.LatLonSphericalBase, bases.LatLonHeightBase, bases.Base
    # test 262 sphericalTrigonometry: sphericalTrigonometry.LatLon, sphericalBase.LatLonSphericalBase, bases.LatLonHeightBase, bases.Base

    # all tests.py tests passed (Python 3.6.0 64bit)
