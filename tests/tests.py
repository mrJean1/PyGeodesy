
# -*- coding: utf-8 -*-

# Tests for various geodesy modules.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# and <http://www.movable-type.co.uk/scripts/latlong.html>.

import sys
from os.path import basename, dirname
try:
    import geodesy as _  # PYCHOK expected
except ImportError:
    # extend sys.path to ../.. directory
    sys.path.insert(0, dirname(dirname(__file__)))
from inspect import isclass, isfunction, ismethod, ismodule

from geodesy import R_M, R_NM, Datums, F_D, F_DM, F_DMS, F_RAD, \
                    compassDMS, compassPoint, degrees, fStr, \
                    lonDMS, normDMS, parseDMS, parse3llh, \
                    precision, toDMS

__all__ = ('Tests',)
__version__ = '16.12.07'

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
        v = sys.version.split()[0]
        n = self.failed
        if n:
            p = '' if n == 1 else 's'
            k = self.known or ''
            if k:
                k = ', incl. %s KNOWN' % (k,)
            r = '(%.1f%%) FAILED%s' % (100.0 * n / self.total, k)
        else:
            n, p, r = 'all', 's', 'passed'
        self.printf('%s %s test%s %s (Python %s)', n, self._name, p, r, v, nl=nl)

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

    def testBases(self, LatLon):
        # bases module tests
        p = LatLon(50.06632, -5.71475)
        self.test('lat, lon', p, '50.06632°N, 005.71475°W')
        q = LatLon('50°03′59″N', """005°42'53"W""")
        self.test('lat, lon', q, '50.066389°N, 005.714722°W')

        p = LatLon(52.205, 0.119)
        q = LatLon(52.205, 0.119)
        self.test('equals', p.equals(q), 'True')

        p = LatLon(51.4778, -0.0016)
        precision(F_DMS, 0)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W''')
        self.test('toStr', p.toStr(F_D), '51.4778°N, 000.0016°W')
        p = LatLon(51.4778, -0.0016, 42)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W, +42.00m''')

    def testDatum(self, geodesy):
        # datum module tests
        E = geodesy.Ellipsoid(1000, 1000, 0, name='TestEllipsiod')
        self.test('ellipsoid', E is geodesy.Ellipsoids.TestEllipsiod, 'True')
#       print(Ellipsoid())

        T = geodesy.Transform(name='TestTransform')
        self.test('transform', T is geodesy.Transforms.TestTransform, 'True')
#       print(Transform())

        D = geodesy.Datum(E, T, name='TestDatum')
        self.test('datum', D is Datums.TestDatum, 'True')
#       print(Datum())

        R, fmt = geodesy.Ellipsoids.WGS84.R, '%.5f'
        self.test('meanR', R, fmt % (R_M,), fmt=fmt)

        E = geodesy.Datums.WGS84.ellipsoid
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = (E.toStr(prec=10),
            'A=%r, e=%.10f, f=1/%.10f, n=%.10f(%.10e)' % (E.A, E.e, 1/E.f, E.n, e),
            'Alpha6=%r' % (E.Alpha6,),
            'Beta6=%r' % (E.Beta6,))
        self.test('WGS84', t[0], "a=6378137.0, b=6356752.3142499998, f=0.0033528107, e2=0.00669438, e22=0.0067394967, R=6371008.7714166669, Rm=6367435.6797186071, name='WGS84'")
        self.test('WGS84', t[1], "A=6367449.145823415, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13)")
        self.test('WGS84', t[2], "Alpha6=(0, 0.0008377318206244698, 7.608527773572307e-07, 1.1976455033294527e-09, 2.4291706072013587e-12, 5.711757677865804e-15, 1.4911177312583895e-17)")
        self.test('WGS84', t[3], "Beta6=(0, 0.0008377321640579486, 5.905870152220203e-08, 1.6734826652839968e-10, 2.1647980400627059e-13, 3.7879780461686053e-16, 7.2487488906941545e-19)")

    def testDMS(self):
        # dms module tests
        self.test('parseDMS', parseDMS(  '0.0°'), '0.0')
        self.test('parseDMS', parseDMS(    '0°'), '0.0')
        self.test('parseDMS', parseDMS('''000°00'00"'''),   '0.0')
        self.test('parseDMS', parseDMS('''000°00'00.0"'''), '0.0')
        self.test('parseDMS', parseDMS('''000° 00'00"'''),    '0.0')
        self.test('parseDMS', parseDMS('''000°00 ' 00.0"'''), '0.0')

        x = parse3llh('000° 00′ 05.31″W, 51° 28′ 40.12″ N')
        x = ', '.join('%.6f' % a for a in x)  # XXX fStr
        self.test('parse3llh', x, '51.477811, -0.001475, 0.000000')

        for a, x in (((),            '''45°45'45.36"'''),
                     ((F_D, None),     '45.7626°'),
                     ((F_DM, None),    "45°45.756'"),
                     ((F_DMS, None), '''45°45'45.36"'''),
                     ((F_D, 6),     '45.7626°'),
                     ((F_DM, -4),   "45°45.7560'"),
                     ((F_DMS, 2), '''45°45'45.36"''')):
            self.test('toDMS', toDMS(45.76260, *a), x)

        for a, x in (((1,),   'N'),
                     ((0,),   'N'),
                     ((-1,),  'N'),
                     ((359,), 'N'),
                     ((24,),   'NNE'),
                     ((24, 1), 'N'),
                     ((24, 2), 'NE'),
                     ((24, 3), 'NNE'),
                     ((226,),   'SW'),
                     ((226, 1), 'W'),
                     ((226, 2), 'SW'),
                     ((226, 3), 'SW'),
                     ((237,),   'WSW'),
                     ((237, 1), 'W'),
                     ((237, 2), 'SW'),
                     ((237, 3), 'WSW')):
            self.test('compassPoint', compassPoint(*a), x)

    def testEllipsoidal(self, LatLon, Nvector=None, Cartesian=None):
        # ellipsoidal modules tests
        p = LatLon(51.4778, -0.0016, 0, Datums.WGS84)
        d = p.convertDatum(Datums.OSGB36)
        self.test('convertDatum', d, '51.477284°N, 000.00002°E, -45.91m')  # 51.4773°N, 000.0000°E, -45.91m
        self.test('convertDatum', d.toStr(F_D, prec=4), '51.4773°N, 000.0°E, -45.91m')

        if Cartesian:
            c = Cartesian(3980581, 97, 4966825)
            n = c.toNvector()  # {x: 0.6228, y: 0.0000, z: 0.7824, h: 0.0000}  # XXX height
            self.test('toNVector', n.toStr(4), '(0.6228, 0.0, 0.7824, +0.24)')
            c = n.toCartesian()
            self.test('toCartesian', c.toStr(0), '[3980581, 97, 4966825]')

        if Nvector:
            n = Nvector(0.5, 0.5, 0.7071)
            c = n.toCartesian()  # [3194434, 3194434, 4487327]
            self.test('toCartesian', c, '[3194434.411, 3194434.411, 4487326.82]')
            p = c.toLatLon()  # 45.0°N, 45.0°E
            self.test('toLatLon', p.toStr('d', 2), '45.0°N, 045.0°E, +0.00m')  # 45.0°N, 45.0°E

            self.test('Nvector', n, '(0.5, 0.5, 0.7071)')
            n = Nvector(0.5, 0.5, 0.7071, 1).toStr(3)
            self.test('Nvector', n, '(0.5, 0.5, 0.707, +1.00)')

    def testLatLon(self, LatLon):
        # basic LatLon class tests
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
            for a in ('lat', 'lon', 'height', 'h'):
                if hasattr(ll, a):
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

    def testSpherical(self, LatLon, Nvector=None):
        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        i = p.intermediateTo(q, 0.25)
        self.test('intermediateTo', i, '51.372084°N, 000.707337°E')

        if hasattr(LatLon, 'intermediateChordTo'):
            i = p.intermediateChordTo(q, 0.25)
            self.test('intermediateChordTo', i, '51.372294°N, 000.707192°E')

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
        self.test('VincentyError' + n, t, 'LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)')

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

    from geodesy import datum, dms, mgrs, osgr, utils, \
                        ellipsoidalNvector, ellipsoidalVincenty, \
                        sphericalNvector, sphericalTrigonometry, \
                        nvector, vector3d, utm
    import geodesy

    t = Tests(__file__, __version__)
    # check that __all__ names exist in each module
    t.testModule(geodesy, 'geodesy')
    for m in (datum, dms, mgrs, osgr, utils,
              ellipsoidalNvector, ellipsoidalVincenty,
              sphericalNvector, sphericalTrigonometry,
              nvector, vector3d, utm):
        t.testModule(m)
    t.testLatLonAttr(ellipsoidalNvector, ellipsoidalVincenty,
                     sphericalNvector, sphericalTrigonometry)
    t.results(nl=1)
    t.exit()

    # Typical test results (on MacOS X)

    # testing tests.py version 16.12.07

    # testing __init__.py version 16.11.11
    # test 1 geodesy.Datum() class (geodesy.datum): True
    # test 2 geodesy.Datums attribute (geodesy.datum): True
    # test 3 geodesy.EPS float: True
    # test 4 geodesy.EPS1 float: True
    # test 5 geodesy.EPS2 float: True
    # test 6 geodesy.Ellipsoid() class (geodesy.datum): True
    # test 7 geodesy.Ellipsoids attribute (geodesy.datum): True
    # test 8 geodesy.F_D str: True
    # test 9 geodesy.F_DM str: True
    # test 10 geodesy.F_DMS str: True
    # test 11 geodesy.F_RAD str: True
    # test 12 geodesy.Mgrs() class (geodesy.mgrs): True
    # test 13 geodesy.Osgr() class (geodesy.osgr): True
    # test 14 geodesy.PI float: True
    # test 15 geodesy.PI2 float: True
    # test 16 geodesy.PI_2 float: True
    # test 17 geodesy.R_KM float: True
    # test 18 geodesy.R_M float: True
    # test 19 geodesy.R_NM float: True
    # test 20 geodesy.R_SM float: True
    # test 21 geodesy.S_DEG str: True
    # test 22 geodesy.S_MIN str: True
    # test 23 geodesy.S_SEC str: True
    # test 24 geodesy.S_SEP str: True
    # test 25 geodesy.Transform() class (geodesy.datum): True
    # test 26 geodesy.Transforms attribute (geodesy.datum): True
    # test 27 geodesy.Utm() class (geodesy.utm): True
    # test 28 geodesy.VincentyError() class (geodesy.ellipsoidalVincenty): True
    # test 29 geodesy.bearingDMS() function (geodesy.dms): True
    # test 30 geodesy.cbrt() function (geodesy.utils): True
    # test 31 geodesy.compassDMS() function (geodesy.dms): True
    # test 32 geodesy.compassPoint() function (geodesy.dms): True
    # test 33 geodesy.degrees attribute (math): True
    # test 34 geodesy.degrees180() function (geodesy.utils): True
    # test 35 geodesy.degrees360() function (geodesy.utils): True
    # test 36 geodesy.degrees90() function (geodesy.utils): True
    # test 37 geodesy.ellipsoidalNvector module: True
    # test 38 geodesy.ellipsoidalVincenty module: True
    # test 39 geodesy.fStr() function (geodesy.utils): True
    # test 40 geodesy.fdot() function (geodesy.utils): True
    # test 41 geodesy.fdot3() function (geodesy.utils): True
    # test 42 geodesy.fsum attribute (math): True
    # test 43 geodesy.halfs() function (geodesy.utils): True
    # test 44 geodesy.hypot1() function (geodesy.utils): True
    # test 45 geodesy.hypot3() function (geodesy.utils): True
    # test 46 geodesy.isint() function (geodesy.utils): True
    # test 47 geodesy.isscalar() function (geodesy.utils): True
    # test 48 geodesy.latDMS() function (geodesy.dms): True
    # test 49 geodesy.len2() function (geodesy.utils): True
    # test 50 geodesy.lonDMS() function (geodesy.dms): True
    # test 51 geodesy.map2() function (geodesy.utils): True
    # test 52 geodesy.normDMS() function (geodesy.dms): True
    # test 53 geodesy.parse3llh() function (geodesy.dms): True
    # test 54 geodesy.parseDMS() function (geodesy.dms): True
    # test 55 geodesy.parseMGRS() function (geodesy.mgrs): True
    # test 56 geodesy.parseOSGR() function (geodesy.osgr): True
    # test 57 geodesy.parseUTM() function (geodesy.utm): True
    # test 58 geodesy.precision() function (geodesy.dms): True
    # test 59 geodesy.radians attribute (math): True
    # test 60 geodesy.radiansPI() function (geodesy.utils): True
    # test 61 geodesy.radiansPI_2() function (geodesy.utils): True
    # test 62 geodesy.sin_2() function (geodesy.utils): True
    # test 63 geodesy.sphericalNvector module: True
    # test 64 geodesy.sphericalTrigonometry module: True
    # test 65 geodesy.tanPI_2_2() function (geodesy.utils): True
    # test 66 geodesy.toDMS() function (geodesy.dms): True
    # test 67 geodesy.toMgrs() function (geodesy.mgrs): True
    # test 68 geodesy.toOsgr() function (geodesy.osgr): True
    # test 69 geodesy.toUtm() function (geodesy.utm): True
    # test 70 geodesy.wrap180() function (geodesy.utils): True
    # test 71 geodesy.wrap90() function (geodesy.utils): True
    # test 72 geodesy.wrapPI() function (geodesy.utils): True
    # test 73 geodesy.wrapPI2() function (geodesy.utils): True
    # test 74 geodesy.wrapPI_2() function (geodesy.utils): True

    # testing datum.py version 16.12.06
    # test 75 datum.Datum() class: True
    # test 76 datum.Datums attribute: True
    # test 77 datum.Ellipsoid() class: True
    # test 78 datum.Ellipsoids attribute: True
    # test 79 datum.R_KM float: True
    # test 80 datum.R_M float: True
    # test 81 datum.R_NM float: True
    # test 82 datum.R_SM float: True
    # test 83 datum.Transform() class: True
    # test 84 datum.Transforms attribute: True

    # testing dms.py version 16.10.10
    # test 85 dms.F_D str: True
    # test 86 dms.F_DM str: True
    # test 87 dms.F_DMS str: True
    # test 88 dms.F_RAD str: True
    # test 89 dms.S_DEG str: True
    # test 90 dms.S_MIN str: True
    # test 91 dms.S_SEC str: True
    # test 92 dms.S_SEP str: True
    # test 93 dms.bearingDMS() function: True
    # test 94 dms.compassDMS() function: True
    # test 95 dms.compassPoint() function: True
    # test 96 dms.latDMS() function: True
    # test 97 dms.lonDMS() function: True
    # test 98 dms.normDMS() function: True
    # test 99 dms.parse3llh() function: True
    # test 100 dms.parseDMS() function: True
    # test 101 dms.precision() function: True
    # test 102 dms.toDMS() function: True

    # testing mgrs.py version 16.11.11
    # test 103 mgrs.Mgrs() class: True
    # test 104 mgrs.parseMGRS() function: True
    # test 105 mgrs.toMgrs() function: True

    # testing osgr.py version 16.11.11
    # test 106 osgr.Osgr() class: True
    # test 107 osgr.parseOSGR() function: True
    # test 108 osgr.toOsgr() function: True

    # testing utils.py version 16.12.04
    # test 109 utils.EPS float: True
    # test 110 utils.EPS1 float: True
    # test 111 utils.EPS2 float: True
    # test 112 utils.PI float: True
    # test 113 utils.PI2 float: True
    # test 114 utils.PI_2 float: True
    # test 115 utils.cbrt() function: True
    # test 116 utils.degrees attribute (math): True
    # test 117 utils.degrees180() function: True
    # test 118 utils.degrees360() function: True
    # test 119 utils.degrees90() function: True
    # test 120 utils.fStr() function: True
    # test 121 utils.fdot() function: True
    # test 122 utils.fdot3() function: True
    # test 123 utils.fsum attribute (math): True
    # test 124 utils.halfs() function: True
    # test 125 utils.hypot1() function: True
    # test 126 utils.hypot3() function: True
    # test 127 utils.isint() function: True
    # test 128 utils.isscalar() function: True
    # test 129 utils.len2() function: True
    # test 130 utils.map2() function: True
    # test 131 utils.radians attribute (math): True
    # test 132 utils.radiansPI() function: True
    # test 133 utils.radiansPI_2() function: True
    # test 134 utils.sin_2() function: True
    # test 135 utils.tanPI_2_2() function: True
    # test 136 utils.wrap180() function: True
    # test 137 utils.wrap90() function: True
    # test 138 utils.wrapPI() function: True
    # test 139 utils.wrapPI2() function: True
    # test 140 utils.wrapPI_2() function: True

    # testing ellipsoidalNvector.py version 16.11.11
    # test 141 ellipsoidalNvector.Cartesian() class: True
    # test 142 ellipsoidalNvector.LatLon() class: True
    # test 143 ellipsoidalNvector.Ned() class: True
    # test 144 ellipsoidalNvector.Nvector() class: True
    # test 145 ellipsoidalNvector.meanOf() function: True
    # test 146 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.py version 16.11.28
    # test 147 ellipsoidalVincenty.Cartesian() class: True
    # test 148 ellipsoidalVincenty.LatLon() class: True
    # test 149 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.py version 16.11.20
    # test 150 sphericalNvector.LatLon() class: True
    # test 151 sphericalNvector.areaOf() function: True
    # test 152 sphericalNvector.intersection() function: True
    # test 153 sphericalNvector.meanOf() function: True
    # test 154 sphericalNvector.triangulate() function: True
    # test 155 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.py version 16.11.11
    # test 156 sphericalTrigonometry.LatLon() class: True
    # test 157 sphericalTrigonometry.meanOf() function: True

    # testing nvector.py version 16.11.11
    # test 158 nvector.NorthPole attribute: True
    # test 159 nvector.Nvector() class: True
    # test 160 nvector.SouthPole attribute: True
    # test 161 nvector.sumOf() function: True

    # testing vector3d.py version 16.11.30
    # test 162 vector3d.Vector3d() class: True
    # test 163 vector3d.sumOf() function: True

    # testing utm.py version 16.12.06
    # test 164 utm.Utm() class: True
    # test 165 utm.parseUTM() function: True
    # test 166 utm.toUtm() function: True

    # testing LatLon.attrs version 16.12.07
    # test 167 Top() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 168 _Nv attribute: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 169 _alter() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 170 _datum attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 171 _direct() method: geodesy.ellipsoidalVincenty
    # test 172 _epsilon float: geodesy.ellipsoidalVincenty
    # test 173 _gc3() method: geodesy.sphericalNvector
    # test 174 _height int: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 175 _inverse() method: geodesy.ellipsoidalVincenty
    # test 176 _iterations int: geodesy.ellipsoidalVincenty
    # test 177 _lat attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 178 _lon attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 179 _osgr attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 180 _r3 attribute: geodesy.ellipsoidalNvector
    # test 181 _rhumb3() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 182 _rotation3() method: geodesy.ellipsoidalNvector
    # test 183 _update() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 184 _utm attribute: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 185 _v3d attribute: geodesy.sphericalTrigonometry
    # test 186 bearingTo() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 187 convertDatum() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 188 copy() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 189 crossTrackDistanceTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 190 crossingParallels() method: geodesy.sphericalTrigonometry
    # test 191 datum property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 192 deltaTo() method: geodesy.ellipsoidalNvector
    # test 193 destination() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 194 destination2() method: geodesy.ellipsoidalVincenty
    # test 195 destinationNed() method: geodesy.ellipsoidalNvector
    # test 196 destinationPoint() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 197 distanceTo() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 198 distanceTo3() method: geodesy.ellipsoidalVincenty
    # test 199 ellipsoid() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 200 ellipsoids() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 201 enclosedBy() method: geodesy.sphericalNvector
    # test 202 epsilon property: geodesy.ellipsoidalVincenty
    # test 203 equals() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 204 finalBearingOn() method: geodesy.ellipsoidalVincenty
    # test 205 finalBearingTo() method: geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 206 greatCircle() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 207 greatCircleTo() method: geodesy.sphericalNvector
    # test 208 height property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry, geodesy.sphericalTrigonometry
    # test 209 initialBearingTo() method: geodesy.ellipsoidalVincenty
    # test 210 intermediateChordTo() method: geodesy.sphericalNvector
    # test 211 intermediatePointDirectlyTo() method: geodesy.sphericalNvector
    # test 212 intermediatePointTo() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 213 intermediateTo() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 214 intersection() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 215 isEnclosedBy() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 216 isWithin() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 217 isWithinExtent() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 218 iterations property: geodesy.ellipsoidalVincenty
    # test 219 lat property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry, geodesy.sphericalTrigonometry
    # test 220 lon property: geodesy.ellipsoidalNvector, geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalNvector, geodesy.sphericalTrigonometry, geodesy.sphericalTrigonometry
    # test 221 maxLat() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 222 maxLatitude() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 223 midpointTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 224 minLat() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 225 minLatitude() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 226 nearestOn() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 227 nearestPointOnSegment() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 228 notImplemented() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 229 others() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 230 parse() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 231 rhumbBearingTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 232 rhumbDistanceTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 233 rhumbMidpointTo() method: geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 234 to3xyz() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 235 to4xyzh() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 236 toCartesian() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 237 toDatum() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 238 toNvector() method: geodesy.ellipsoidalNvector, geodesy.sphericalNvector
    # test 239 toOsgr() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 240 toStr() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 241 toStr2() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 242 toUtm() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty
    # test 243 toVector3d() method: geodesy.sphericalTrigonometry
    # test 244 toradians() method: geodesy.ellipsoidalNvector, geodesy.ellipsoidalVincenty, geodesy.sphericalNvector, geodesy.sphericalTrigonometry
    # test 245 triangulate() method: geodesy.sphericalNvector
    # test 246 trilaterate() method: geodesy.sphericalNvector

    # testing LatLon.mro version 16.12.07
    # test 247 geodesy.ellipsoidalNvector: geodesy.ellipsoidalNvector.LatLon, geodesy.nvector._LatLonNvectorBase, geodesy.ellipsoidalBase._LatLonHeightDatumBase, geodesy.bases._LatLonHeightBase, geodesy.bases._Base
    # test 248 geodesy.ellipsoidalVincenty: geodesy.ellipsoidalVincenty.LatLon, geodesy.ellipsoidalBase._LatLonHeightDatumBase, geodesy.bases._LatLonHeightBase, geodesy.bases._Base
    # test 249 geodesy.sphericalNvector: geodesy.sphericalNvector.LatLon, geodesy.nvector._LatLonNvectorBase, geodesy.sphericalBase._LatLonSphericalBase, geodesy.bases._LatLonHeightBase, geodesy.bases._Base
    # test 250 geodesy.sphericalTrigonometry: geodesy.sphericalTrigonometry.LatLon, geodesy.sphericalBase._LatLonSphericalBase, geodesy.bases._LatLonHeightBase, geodesy.bases._Base

    # all tests.py tests passed (Python 2.7.10)

    # testing tests.py version 16.12.07

    # testing __init__.py version 16.11.11
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
    # test 12 geodesy.Mgrs() class (mgrs): True
    # test 13 geodesy.Osgr() class (osgr): True
    # test 14 geodesy.PI float: True
    # test 15 geodesy.PI2 float: True
    # test 16 geodesy.PI_2 float: True
    # test 17 geodesy.R_KM float: True
    # test 18 geodesy.R_M float: True
    # test 19 geodesy.R_NM float: True
    # test 20 geodesy.R_SM float: True
    # test 21 geodesy.S_DEG str: True
    # test 22 geodesy.S_MIN str: True
    # test 23 geodesy.S_SEC str: True
    # test 24 geodesy.S_SEP str: True
    # test 25 geodesy.Transform() class (datum): True
    # test 26 geodesy.Transforms attribute (datum): True
    # test 27 geodesy.Utm() class (utm): True
    # test 28 geodesy.VincentyError() class (ellipsoidalVincenty): True
    # test 29 geodesy.bearingDMS() function (dms): True
    # test 30 geodesy.cbrt() function (utils): True
    # test 31 geodesy.compassDMS() function (dms): True
    # test 32 geodesy.compassPoint() function (dms): True
    # test 33 geodesy.degrees attribute (math): True
    # test 34 geodesy.degrees180() function (utils): True
    # test 35 geodesy.degrees360() function (utils): True
    # test 36 geodesy.degrees90() function (utils): True
    # test 37 geodesy.ellipsoidalNvector module: True
    # test 38 geodesy.ellipsoidalVincenty module: True
    # test 39 geodesy.fStr() function (utils): True
    # test 40 geodesy.fdot() function (utils): True
    # test 41 geodesy.fdot3() function (utils): True
    # test 42 geodesy.fsum attribute (math): True
    # test 43 geodesy.halfs() function (utils): True
    # test 44 geodesy.hypot1() function (utils): True
    # test 45 geodesy.hypot3() function (utils): True
    # test 46 geodesy.isint() function (utils): True
    # test 47 geodesy.isscalar() function (utils): True
    # test 48 geodesy.latDMS() function (dms): True
    # test 49 geodesy.len2() function (utils): True
    # test 50 geodesy.lonDMS() function (dms): True
    # test 51 geodesy.map2() function (utils): True
    # test 52 geodesy.normDMS() function (dms): True
    # test 53 geodesy.parse3llh() function (dms): True
    # test 54 geodesy.parseDMS() function (dms): True
    # test 55 geodesy.parseMGRS() function (mgrs): True
    # test 56 geodesy.parseOSGR() function (osgr): True
    # test 57 geodesy.parseUTM() function (utm): True
    # test 58 geodesy.precision() function (dms): True
    # test 59 geodesy.radians attribute (math): True
    # test 60 geodesy.radiansPI() function (utils): True
    # test 61 geodesy.radiansPI_2() function (utils): True
    # test 62 geodesy.sin_2() function (utils): True
    # test 63 geodesy.sphericalNvector module: True
    # test 64 geodesy.sphericalTrigonometry module: True
    # test 65 geodesy.tanPI_2_2() function (utils): True
    # test 66 geodesy.toDMS() function (dms): True
    # test 67 geodesy.toMgrs() function (mgrs): True
    # test 68 geodesy.toOsgr() function (osgr): True
    # test 69 geodesy.toUtm() function (utm): True
    # test 70 geodesy.wrap180() function (utils): True
    # test 71 geodesy.wrap90() function (utils): True
    # test 72 geodesy.wrapPI() function (utils): True
    # test 73 geodesy.wrapPI2() function (utils): True
    # test 74 geodesy.wrapPI_2() function (utils): True

    # testing datum.py version 16.12.06
    # test 75 datum.Datum() class: True
    # test 76 datum.Datums attribute: True
    # test 77 datum.Ellipsoid() class: True
    # test 78 datum.Ellipsoids attribute: True
    # test 79 datum.R_KM float: True
    # test 80 datum.R_M float: True
    # test 81 datum.R_NM float: True
    # test 82 datum.R_SM float: True
    # test 83 datum.Transform() class: True
    # test 84 datum.Transforms attribute: True

    # testing dms.py version 16.10.10
    # test 85 dms.F_D str: True
    # test 86 dms.F_DM str: True
    # test 87 dms.F_DMS str: True
    # test 88 dms.F_RAD str: True
    # test 89 dms.S_DEG str: True
    # test 90 dms.S_MIN str: True
    # test 91 dms.S_SEC str: True
    # test 92 dms.S_SEP str: True
    # test 93 dms.bearingDMS() function: True
    # test 94 dms.compassDMS() function: True
    # test 95 dms.compassPoint() function: True
    # test 96 dms.latDMS() function: True
    # test 97 dms.lonDMS() function: True
    # test 98 dms.normDMS() function: True
    # test 99 dms.parse3llh() function: True
    # test 100 dms.parseDMS() function: True
    # test 101 dms.precision() function: True
    # test 102 dms.toDMS() function: True

    # testing mgrs.py version 16.11.11
    # test 103 mgrs.Mgrs() class: True
    # test 104 mgrs.parseMGRS() function: True
    # test 105 mgrs.toMgrs() function: True

    # testing osgr.py version 16.11.11
    # test 106 osgr.Osgr() class: True
    # test 107 osgr.parseOSGR() function: True
    # test 108 osgr.toOsgr() function: True

    # testing utils.py version 16.12.04
    # test 109 utils.EPS float: True
    # test 110 utils.EPS1 float: True
    # test 111 utils.EPS2 float: True
    # test 112 utils.PI float: True
    # test 113 utils.PI2 float: True
    # test 114 utils.PI_2 float: True
    # test 115 utils.cbrt() function: True
    # test 116 utils.degrees attribute (math): True
    # test 117 utils.degrees180() function: True
    # test 118 utils.degrees360() function: True
    # test 119 utils.degrees90() function: True
    # test 120 utils.fStr() function: True
    # test 121 utils.fdot() function: True
    # test 122 utils.fdot3() function: True
    # test 123 utils.fsum attribute (math): True
    # test 124 utils.halfs() function: True
    # test 125 utils.hypot1() function: True
    # test 126 utils.hypot3() function: True
    # test 127 utils.isint() function: True
    # test 128 utils.isscalar() function: True
    # test 129 utils.len2() function: True
    # test 130 utils.map2() function: True
    # test 131 utils.radians attribute (math): True
    # test 132 utils.radiansPI() function: True
    # test 133 utils.radiansPI_2() function: True
    # test 134 utils.sin_2() function: True
    # test 135 utils.tanPI_2_2() function: True
    # test 136 utils.wrap180() function: True
    # test 137 utils.wrap90() function: True
    # test 138 utils.wrapPI() function: True
    # test 139 utils.wrapPI2() function: True
    # test 140 utils.wrapPI_2() function: True

    # testing ellipsoidalNvector.py version 16.11.11
    # test 141 ellipsoidalNvector.Cartesian() class: True
    # test 142 ellipsoidalNvector.LatLon() class: True
    # test 143 ellipsoidalNvector.Ned() class: True
    # test 144 ellipsoidalNvector.Nvector() class: True
    # test 145 ellipsoidalNvector.meanOf() function: True
    # test 146 ellipsoidalNvector.toNed() function: True

    # testing ellipsoidalVincenty.py version 16.11.28
    # test 147 ellipsoidalVincenty.Cartesian() class: True
    # test 148 ellipsoidalVincenty.LatLon() class: True
    # test 149 ellipsoidalVincenty.VincentyError() class: True

    # testing sphericalNvector.py version 16.11.20
    # test 150 sphericalNvector.LatLon() class: True
    # test 151 sphericalNvector.areaOf() function: True
    # test 152 sphericalNvector.intersection() function: True
    # test 153 sphericalNvector.meanOf() function: True
    # test 154 sphericalNvector.triangulate() function: True
    # test 155 sphericalNvector.trilaterate() function: True

    # testing sphericalTrigonometry.py version 16.11.11
    # test 156 sphericalTrigonometry.LatLon() class: True
    # test 157 sphericalTrigonometry.meanOf() function: True

    # testing nvector.py version 16.11.11
    # test 158 nvector.NorthPole attribute: True
    # test 159 nvector.Nvector() class: True
    # test 160 nvector.SouthPole attribute: True
    # test 161 nvector.sumOf() function: True

    # testing vector3d.py version 16.11.30
    # test 162 vector3d.Vector3d() class: True
    # test 163 vector3d.sumOf() function: True

    # testing utm.py version 16.12.06
    # test 164 utm.Utm() class: True
    # test 165 utm.parseUTM() function: True
    # test 166 utm.toUtm() function: True

    # testing LatLon.attrs version 16.12.07
    # test 167 Top() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 168 _Nv attribute: ellipsoidalNvector, sphericalNvector
    # test 169 _alter() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 170 _datum attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 171 _direct() function: ellipsoidalVincenty
    # test 172 _epsilon float: ellipsoidalVincenty
    # test 173 _gc3() function: sphericalNvector
    # test 174 _height int: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 175 _inverse() function: ellipsoidalVincenty
    # test 176 _iterations int: ellipsoidalVincenty
    # test 177 _lat attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 178 _lon attribute: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 179 _osgr attribute: ellipsoidalNvector, ellipsoidalVincenty
    # test 180 _r3 attribute: ellipsoidalNvector
    # test 181 _rhumb3() function: sphericalNvector, sphericalTrigonometry
    # test 182 _rotation3() function: ellipsoidalNvector
    # test 183 _update() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 184 _utm attribute: ellipsoidalNvector, ellipsoidalVincenty
    # test 185 _v3d attribute: sphericalTrigonometry
    # test 186 bearingTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 187 convertDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 188 copy() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 189 crossTrackDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 190 crossingParallels() function: sphericalTrigonometry
    # test 191 datum property: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 192 deltaTo() function: ellipsoidalNvector
    # test 193 destination() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 194 destination2() function: ellipsoidalVincenty
    # test 195 destinationNed() function: ellipsoidalNvector
    # test 196 destinationPoint() function: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 197 distanceTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 198 distanceTo3() function: ellipsoidalVincenty
    # test 199 ellipsoid() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 200 ellipsoids() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 201 enclosedBy() function: sphericalNvector
    # test 202 epsilon property: ellipsoidalVincenty
    # test 203 equals() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 204 finalBearingOn() function: ellipsoidalVincenty
    # test 205 finalBearingTo() function: ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 206 greatCircle() function: sphericalNvector, sphericalTrigonometry
    # test 207 greatCircleTo() function: sphericalNvector
    # test 208 height property: ellipsoidalNvector, ellipsoidalNvector, ellipsoidalVincenty, ellipsoidalVincenty, sphericalNvector, sphericalNvector, sphericalTrigonometry, sphericalTrigonometry
    # test 209 initialBearingTo() function: ellipsoidalVincenty
    # test 210 intermediateChordTo() function: sphericalNvector
    # test 211 intermediatePointDirectlyTo() function: sphericalNvector
    # test 212 intermediatePointTo() function: ellipsoidalNvector, sphericalNvector
    # test 213 intermediateTo() function: ellipsoidalNvector, sphericalNvector, sphericalTrigonometry
    # test 214 intersection() function: sphericalNvector, sphericalTrigonometry
    # test 215 isEnclosedBy() function: sphericalNvector, sphericalTrigonometry
    # test 216 isWithin() function: sphericalNvector, sphericalTrigonometry
    # test 217 isWithinExtent() function: sphericalNvector, sphericalTrigonometry
    # test 218 iterations property: ellipsoidalVincenty
    # test 219 lat property: ellipsoidalNvector, ellipsoidalNvector, ellipsoidalVincenty, ellipsoidalVincenty, sphericalNvector, sphericalNvector, sphericalTrigonometry, sphericalTrigonometry
    # test 220 lon property: ellipsoidalNvector, ellipsoidalNvector, ellipsoidalVincenty, ellipsoidalVincenty, sphericalNvector, sphericalNvector, sphericalTrigonometry, sphericalTrigonometry
    # test 221 maxLat() function: sphericalNvector, sphericalTrigonometry
    # test 222 maxLatitude() function: sphericalNvector, sphericalTrigonometry
    # test 223 midpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 224 minLat() function: sphericalNvector, sphericalTrigonometry
    # test 225 minLatitude() function: sphericalNvector, sphericalTrigonometry
    # test 226 nearestOn() function: sphericalNvector, sphericalTrigonometry
    # test 227 nearestPointOnSegment() function: sphericalNvector, sphericalTrigonometry
    # test 228 notImplemented() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 229 others() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 230 parse() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 231 rhumbBearingTo() function: sphericalNvector, sphericalTrigonometry
    # test 232 rhumbDistanceTo() function: sphericalNvector, sphericalTrigonometry
    # test 233 rhumbMidpointTo() function: sphericalNvector, sphericalTrigonometry
    # test 234 to3xyz() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 235 to4xyzh() function: ellipsoidalNvector, sphericalNvector
    # test 236 toCartesian() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 237 toDatum() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 238 toNvector() function: ellipsoidalNvector, sphericalNvector
    # test 239 toOsgr() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 240 toStr() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 241 toStr2() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 242 toUtm() function: ellipsoidalNvector, ellipsoidalVincenty
    # test 243 toVector3d() function: sphericalTrigonometry
    # test 244 toradians() function: ellipsoidalNvector, ellipsoidalVincenty, sphericalNvector, sphericalTrigonometry
    # test 245 triangulate() function: sphericalNvector
    # test 246 trilaterate() function: sphericalNvector

    # testing LatLon.mro version 16.12.07
    # test 247 ellipsoidalNvector: ellipsoidalNvector.LatLon, nvector._LatLonNvectorBase, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 248 ellipsoidalVincenty: ellipsoidalVincenty.LatLon, ellipsoidalBase._LatLonHeightDatumBase, bases._LatLonHeightBase, bases._Base
    # test 249 sphericalNvector: sphericalNvector.LatLon, nvector._LatLonNvectorBase, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base
    # test 250 sphericalTrigonometry: sphericalTrigonometry.LatLon, sphericalBase._LatLonSphericalBase, bases._LatLonHeightBase, bases._Base

    # all tests.py tests passed (Python 3.5.2)
