
# -*- coding: utf-8 -*-

# Tests for various pygeodesy modules.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# and <http://www.movable-type.co.uk/scripts/latlong.html>.

from os.path import basename, dirname
import sys
try:
    import pygeodesy as _  # PYCHOK expected
except ImportError:
    # extend sys.path to ../.. directory
    sys.path.insert(0, dirname(dirname(__file__)))
from pygeodesy import R_NM, F_D, F_DM, F_DMS, F_RAD, \
                      version as geodesy_version, \
                      degrees, isclockwise, m2NM, normDMS  # PYCHOK expected

from inspect import isclass, isfunction, ismethod, ismodule
from platform import architecture
from time import time

__all__ = ('versions', 'Tests',
           'secs2str')
__version__ = '17.05.29'

try:
    _int = int, long
    _str = basestring
except NameError:  # Python 3+
    _int = int
    _str = str

versions = ' '.join(('PyGeodesy', geodesy_version,
                     'Python', sys.version.split()[0], architecture()[0]))


def secs2str(secs):
    unit = ['sec', 'ms', 'us', 'ps']
    while secs < 1 and len(unit) > 1:
        secs *= 1000.0
        unit.pop(0)
    return '%.3f %s' % (secs, unit[0])


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
       and examples in <http://www.edwilliams.org/avform.htm> or
       elsewhere as indicated.
    '''
    _file   = ''
    _name   = ''
    _prefix = '    '
    _tests  = []
    _time   = 0

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
        self._time = time()

    def errors(self):
        return self.failed - self.known  # new failures

    def exit(self, errors=0):
        sys.exit(min(errors + self.errors(), 99))

    def printf(self, fmt, *args, **kwds):  # nl=0
        nl = '\n' * kwds.get('nl', 0)
        print(nl + self._prefix + (fmt % args))

    def results(self, nl=0):
        n = self.failed
        if n:
            p = '' if n == 1 else 's'
            k = self.known or ''
            if k:
                k = ', incl. %s KNOWN' % (k,)
            r = '(%.1f%%) FAILED%s' % (100.0 * n / self.total, k)
        else:
            n, p, r = 'all', 's', 'passed'
        s = time() - self._time
        t = '(%s) %s' % (versions, secs2str(s))
        self.printf('%s %s test%s %s %s', n, self._name, p, r, t, nl=nl)

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
        self.printf('testing module %s version %s', basename(module), version, nl=1)

    def testLatLon(self, LatLon, Sph=True):  # MCCABE expected
        # basic LatLon class tests
        p = LatLon(52.20472, 0.14056)
        self.test('isellipsoidal', p.isellipsoidal, str(not Sph))
        self.test('isspherical', p.isspherical, str(Sph))

        self.test('lat/lonDMS', p, '52.20472°N, 000.14056°E')  # 52.20472°N, 000.14056°E
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 3),  '''52°12.283'N, 000°08.434'E''')
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 4),  '''52°12.2832'N, 000°08.4336'E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 0), '''52°12'17"N, 000°08'26"E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 1), '''52°12'17.0"N, 000°08'26.0"E''')
        self.test('lat/lonDMS F_RAD', p.toStr(F_RAD, 6), '0.911144N, 0.002453E')
        q = LatLon(*map(degrees, p.to2ab()))
        self.test('equals', q.equals(p), 'True')

        # <http://www.edwilliams.org/avform.htm#XTE>
        LAX = LatLon(33.+57./60, -(118.+24./60))
        JFK = LatLon(degrees(0.709186), -degrees(1.287762))
        Rav = m2NM(6366710)  # av earth radius in NM

        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        self.test('equals', p.equals(q), 'False')

        if hasattr(LatLon, 'initialBearingTo'):
            b = p.initialBearingTo(q)
            self.test('initialBearingTo', b, '156.1666' if Sph else '156.1106', '%.4f')  # 156.2
            b = p.finalBearingTo(q)
            self.test('finalBearingTo', b, '157.8904' if Sph else '157.8345', '%.4f')
            b = LAX.initialBearingTo(JFK)
            self.test('initialBearingTo', b, '65.8921' if Sph else '65.9335', '%.4f')  # PYCHOK false?  66

        c = p.copy()
        self.test('copy', p.equals(c), 'True')

        if hasattr(LatLon, 'distanceTo'):
            d = p.distanceTo(q)
            self.test('distanceTo', d, '404279.720589' if Sph else '404607.805988', '%.6f')  # 404300
            d = q.distanceTo(p)
            self.test('distanceTo', d, '404279.720589' if Sph else '404607.805988', '%.6f')  # 404300
            d = LAX.distanceTo(JFK, radius=R_NM) if Sph else LAX.distanceTo(JFK)
            self.test('distanceTo', d, '2145' if Sph else '3981601', '%.0f')  # PYCHOK false?

        if hasattr(LatLon, 'intermediateTo'):
            i = p.intermediateTo(q, 0.25)
            self.test('intermediateTo', i, '51.372084°N, 000.707337°E' if Sph
                                      else '51.372294°N, 000.707192°E')
            self.test('intermediateTo', isinstance(i, LatLon), 'True')

            if hasattr(p, 'distanceTo'):
                d = p.distanceTo(q)
                self.test('intermediateTo', d, '404279.721', '%.3f')  # PYCHOK false?

            i = p.intermediateTo(q, 5)
            self.test('intermediateTo+5', i, '35.160975°N, 008.989542°E' if Sph
                                        else '35.560239°N, 008.833512°E')
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo+5', p.distanceTo(i) / d, '5.000', '%.3f')  # PYCHOK false?

            i = p.intermediateTo(q, -4)
            self.test('intermediateTo-4', i, '64.911647°N, 013.726301°W' if Sph
                                        else '64.570387°N, 013.156352°W')
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo-4', p.distanceTo(i) / d, '4.000', '%.3f')  # PYCHOK false?

        if hasattr(LatLon, 'intermediateChordTo'):
            i = p.intermediateChordTo(q, 0.25)
            self.test('intermediateChordTo', i, '51.372294°N, 000.707192°E')
            self.test('intermediateChordTo', isinstance(i, LatLon), 'True')  # PYCHOK false?

        if hasattr(LatLon, 'midpointTo'):
            m = p.midpointTo(q)
            self.test('midpointTo', m, '50.536327°N, 001.274614°E')  # PYCHOK false?  # 50.5363°N, 001.2746°E

        if hasattr(LatLon, 'destination'):
            p = LatLon(51.4778, -0.0015)
            d = p.destination(7794, 300.7)
            self.test('destination', d, '51.513546°N, 000.098345°W' if Sph
                                   else '51.513526°N, 000.098038°W')  # 51.5135°N, 0.0983°W ???
            self.test('destination', d.toStr(F_DMS, 0), '51°30′49″N, 000°05′54″W' if Sph
                                                   else '51°30′49″N, 000°05′53″W')
            d = LAX.destination(100, 66, radius=R_NM) if Sph else LAX.destination(100, 66)
            self.test('destination', d.toStr(F_DM, prec=0), "34°37′N, 116°33′W" if Sph
                                                       else "33°57′N, 118°24′W")
            self.test('destination', d, '34.613647°N, 116.55116°W' if Sph
                                   else '33.950367°N, 118.399012°W')  # PYCHOK false?

        if hasattr(LatLon, 'alongTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.alongTrackDistanceTo(s, 96)
                self.test('alongTrackDistanceTo', d, '62331.59', '%.2f')  # 62331
            except TypeError as x:
                self.test('alongTrackDistanceTo', x, 'type(end) mismatch: int vs sphericalTrigonometry.LatLon')  # PYCHOK false?
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '62331.58', '%.2f')  # PYCHOK false?

            # <http://www.edwilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.alongTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('alongTrackDistanceTo', d, '99.588', '%.3f')  # NM

            # courtesy of Rimvydas Naktinis
            p = LatLon(53.36366, -1.83883)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '-7702.7', '%.1f')

            p = LatLon(53.35423, -1.60881)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '7587.6', '%.1f')  # PYCHOK false?

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.crossTrackDistanceTo(s, 96)
                self.test('crossTrackDistanceTo', d, '-305.67', '%.2f')  # -305.7
            except TypeError as x:
                self.test('crossTrackDistanceTo', x, 'type(end) mismatch: int vs sphericalTrigonometry.LatLon')  # PYCHOK false?
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', '%.2f')  # PYCHOK false?  # -307.5

            # <http://www.edwilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.crossTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('crossTrackDistanceTo', d, '7.4524', '%.4f')  # PYCHOK false? # XXX 7.4512 NM

        if hasattr(LatLon, 'greatCircle'):
            p = LatLon(53.3206, -1.7297)
            gc = p.greatCircle(96.0)
            self.test('greatCircle', gc, '(-0.79408, 0.12856, 0.59406)')  # PYCHOK false?

        if hasattr(LatLon, 'greatCircleTo'):
            p = LatLon(53.3206, -1.7297)
            q = LatLon(53.1887, 0.1334)
            gc = p.greatCircleTo(q)
            self.test('greatCircleTo', gc, '(-0.79408, 0.12859, 0.59406)')  # PYCHOK false?

        if isclockwise:
            f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            self.test('isclockwise', isclockwise(f), 'False')
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            self.test('isclockwise', isclockwise(t), 'True')
            try:
                self.test('isclockwise', isclockwise(t[:2]), ValueError)
            except ValueError as x:
                self.test('isclockwise', x, 'too few points: 2')  # PYCHOK false?

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

    def testVectorial(self, LatLon, Nvector, sumOf):
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
        self.test('length', v.length, '0.99992449715')
        self.test('length', c.length, '1.0')

        class Nv(Nvector):
            pass
        v = Nvector(52.205, 0.119, 0.0)
        s = sumOf((v, c), Vector=Nv, h=0)
        self.test('sumOf', s, '(52.70504, 0.61904, 0.70705)')
        self.test('sumOf', s.__class__.__name__, 'Nv')
        self.test('length', s.length, '52.7134151513')

        c = v.copy()
        self.test('copy', c.equals(v), 'True')
        self.test('length', v.length, '52.2051356286')
        self.test('length', c.length, '52.2051356286')

        if hasattr(LatLon, 'nearestOn'):
            s1 = LatLon(51.0, 1.0)
            s2 = LatLon(51.0, 2.0)
            s = LatLon(51.0, 1.9)
            p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E
            self.test('nearestOn', p.toStr(F_D, prec=4), '51.0004°N, 001.9°E')
            self.test('nearestOn', isinstance(p, LatLon), 'True')

            d = p.distanceTo(s)  # 42.71 m
            self.test('distanceTo', d, '42.712', fmt='%.3f')
            s = LatLon(51.0, 2.1)
            p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E
            self.test('nearestOn', p.toStr(F_D), '51.0°N, 002.0°E')
            self.test('nearestOn', isinstance(p, LatLon), 'True')

            # courtesy AkimboEG on GitHub
            s1 = LatLon(0, 0)
            s2 = LatLon(0, 1)
            s = LatLon(1, 0)
            p = s.nearestOn(s1, s2)  # 0.0°N, 0.0°E
            self.test('nearestOn', p, '00.0°N, 000.0°E')
            self.test('nearestOn', isinstance(p, LatLon), 'True')

            p = LatLon(10, -140).nearestOn(LatLon(0, 20), LatLon(0, 40))
            self.test('nearestOn', p, '00.0°N, 020.0°E')
            self.test('nearestOn', isinstance(p, LatLon), 'True')

        if hasattr(LatLon, 'triangulate'):
            # courtesy of pvezid  Feb 10, 2017
            p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
            self.test('BasseC', p, '47.3038°N, 002.5721°W')
            s = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
            self.test('BasseH', s, '47.311067°N, 002.528617°W')
            t = p.triangulate(7, s, 295)
            self.test('triangulate', t, '47.323667°N, 002.568501°W')
            self.test('triangulate', isinstance(t, LatLon), 'True')


if __name__ == '__main__':

    from pygeodesy import datum, dms, \
                          ellipsoidalNvector, ellipsoidalVincenty, \
                          lcc, mgrs, nvector, osgr, simplify, \
                          sphericalNvector, sphericalTrigonometry, \
                          vector3d, utm, utils  # PYCHOK expected
    import pygeodesy  # PYCHOK expected

    t = Tests(__file__, __version__)
    # check that __all__ names exist in each module
    t.testModule(pygeodesy, 'pygeodesy')
    for m in (datum, dms,
              ellipsoidalNvector, ellipsoidalVincenty,
              lcc, mgrs, nvector, osgr, simplify,
              sphericalNvector, sphericalTrigonometry,
              vector3d, utm, utils):
        t.testModule(m)
    t.testLatLonAttr(ellipsoidalNvector, ellipsoidalVincenty,
                     sphericalNvector, sphericalTrigonometry)
    t.results(nl=1)
    t.exit()
