
# -*- coding: utf-8 -*-

# Test the simplify functions.

__all__ = ('Tests',)
__version__ = '21.02.11'

from base import geographiclib, isPython37, TestsBase

from pygeodesy import EPS, NN, R_M, R_MA, LatLon_, \
                      LatLon2psxy, Numpy2LatLon, Tuple2LatLon, \
                      areaOf, boundsOf, centroidOf, classname, fstr, \
                      isclockwise, isconvex, perimeterOf, points, quadOf

try:
    if isPython37:
        from collections.abc import Sequence
    else:
        from collections import Sequence
except ImportError:
    Sequence = None

try:
    import numpy
except ImportError:
    numpy = None


class Tests(TestsBase):

    def test2(self, pts, npt, psxy=False):

        clas_ = classname(pts) + '.'

        def _test(name, *args, **kwds):
            # prefix class to test name
            self.test(clas_ + name, *args, **kwds)

        if Sequence:  # check abstact base class conformance
            _test('ABC', isinstance(pts, Sequence), True)

        e = pts.epsilon
        _test('epsilon', e, EPS)
        pts.epsilon = 0
        _test('epsilon', pts.epsilon, 0.0)
        pts.epsilon = e

        n = len(pts) // 6  # 0 < some number < len(pts)
        _test('len',    len(pts), len(npt))
        _test('iter',   len(tuple(iter(pts))), len(npt))
        if hasattr(npt, 'shape'):
            _test('shape',  npt.shape, pts.shape)
        _test('slice1', len(pts[:n]), n)
        _test('slice2', type(pts[1:n:2]), type(pts))
        _test('slice3', pts[1:n][0], pts[1])
        _test('str/repr', str(pts), repr(pts))
        if hasattr(pts, 'subset'):
            _test('subset', type(pts.subset(range(n))), type(npt))  # , nt=1)

        for i, p in ((10, (52.224006, -0.707747)),
                     (20, (52.232688, -0.714608)),
                     (30, (52.234375, -0.714348)),
                     (40, (52.237239, -0.712557)),
                     (50, (52.24023,  -0.709919)),
                     (60, (52.240745, -0.707042))):

            if psxy:
                _test('find LL', pts.find(LatLon_(*p)), i)  # coverage
                p = tuple(reversed(p))  # flip to x, y tuple
                _test('find LL', pts.find(LatLon_(*p)), -1)  # coverage

            _test('count', pts.count(p), 1)
            _test('index', pts.index(p), i)
            _test('rfind', pts.rfind(p), i)
            _test('in', p in pts, True)

            p = tuple(reversed(p))
            _test('count', pts.count(p), 0)
            _test('find', pts.find(p), -1)
            _test('rfind', pts.rfind(p), -1)
            _test('not in', p not in pts, True)

        pts = pts[::6]
        for i, p in enumerate(pts):
            _test('enumerate[%s]' % (i,), p, pts[i])
        _test('enumerate[%s]' % ('*',), i, len(pts) - 1)

        i = len(pts)
        for p in reversed(pts):
            i -= 1
            _test('reversed[%s]' % (i,), p, pts[i])

        for i, p in enumerate(pts):
            _test('findall[%s]' % (i,), tuple(pts.findall(p)), (i,))
        _test('findall[%s]' % ('*',), i, len(pts) - 1)

        _test('isNumpy2',  pts.isNumpy2,  pts.__class__ is Numpy2LatLon)
        _test('isPoints2', pts.isPoints2, pts.__class__ is LatLon2psxy)
        _test('isTuple2',  pts.isTuple2,  pts.__class__ is Tuple2LatLon)

    def test3(self, LatLon):
        self.subtitle(points, LatLon=LatLon)

        p = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        self.test('areaOf', areaOf(p, radius=R_MA), '8.811228e+09', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=6), '45.5, 1.5')
        self.test('perimeterOf', perimeterOf(p, radius=R_MA), '2.673633e+05', fmt='%.6e')
        self.test('isclockwise', isclockwise(p), False)
        self.test('isconvex', isconvex(p), True)

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('areaOf', areaOf(p, radius=R_MA), '7.086883e+09', fmt='%.6e')
        self.test('perimeterOf', perimeterOf(p, radius=R_MA), '2.687460e+05', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=6), '0.333333, 0.333333')
        self.test('isclockwise', isclockwise(p), True)
        self.test('isconvex', isconvex(p), True)

        p = LatLon(0, 1), LatLon(1, 2), LatLon(2, 1), LatLon(1, 0)
        self.test('areaOf', areaOf(p, radius=R_M), '2.827856e+10', fmt='%.6e')
        self.test('perimeterOf', perimeterOf(p, radius=R_M), '4.717039e+05', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=6), '1.0, 1.0')
        self.test('isclockwise', isclockwise(p), False)
        self.test('isconvex', isconvex(p), True)

        p = LatLon(45, -70), LatLon(60, 0), LatLon(20, 110), LatLon(80, 170)  # XXX warped?
        self.test('areaOf', areaOf(p, radius=R_M), '2.747297e+13', fmt='%.6e')  # XXX 1.047657e+12
        self.test('perimeterOf', perimeterOf(p, radius=R_M), '2.332643e+07', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=3), '52.113, 102.123')  # XXX '22.536, -164.928'
        self.test('isclockwise', isclockwise(p), False)  # XXX True
        self.test('isconvex', isconvex(p), False)

        p = LatLon(0, 0), LatLon(0, 3), LatLon(3, 3), LatLon(3, 2), \
            LatLon(1, 2), LatLon(1, 1), LatLon(2, 1), LatLon(2, 0)
        self.test('areaOf', areaOf(p, radius=R_M), '8.482014e+10', fmt='%.6e')
        self.test('perimeterOf', perimeterOf(p, radius=R_M), '1.334104e+06', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=3), '1.167, 1.667')
        self.test('isclockwise', isclockwise(p), False)
        self.test('isconvex', isconvex(p), False)

        p = LatLon(-20, -180), LatLon(5, -160), LatLon(0, -60), LatLon(-60, -160)
        self.test('areaOf', areaOf(p, radius=R_M), '5.151974e+13', fmt='%.6e')
        self.test('perimeterOf', perimeterOf(p, radius=R_M), '2.638608e+07', fmt='%.6e')
        self.test('centroidOf', fstr(centroidOf(p), prec=3), '-19.444, -133.333')
        self.test('isclockwise', isclockwise(p), True)
        self.test('isconvex', isconvex(p), True)

        # <https://GeographicLib.SourceForge.io/scripts/geod-calc.html>
        p = (LatLon(-63.1,  -58), LatLon(-72.9,  -74), LatLon(-71.9, -102),
             LatLon(-74.9, -102), LatLon(-74.3, -131), LatLon(-77.5, -163),
             LatLon(-77.4,  163), LatLon(-71.7,  172), LatLon(-65.9,  140),
             LatLon(-65.7,  113), LatLon(-66.6,   88), LatLon(-66.9,   59),
             LatLon(-69.8,   25), LatLon(-70.0,   -4), LatLon(-71.0,  -14),
             LatLon(-77.3,  -33), LatLon(-77.9,  -46), LatLon(-74.7,  -61))  # on/around south pole!
        self.test('areaOf', areaOf(p, radius=R_M), '4.469277e+12', fmt='%.6e', known=True)
        self.test('perimeterOf', perimeterOf(p, radius=R_M), '1.562029e+07', fmt='%.6e', known=True)
        self.test('centroidOf', fstr(centroidOf(p), prec=3), '-71.443, -69.683')  # XXX '-72.112, 92.032'
        self.test('isclockwise', isclockwise(p), True)  # XXX False
        self.test('isconvex', isconvex(p), False)
        self.test('points2', p[0].points2(p)[0], len(p))

        if LatLon is LatLon_:
            b = boundsOf(p)
            self.test('boundsOf', b, '(-77.9, -163.0, -63.1, 172.0)')
            q = str(quadOf(*b)).replace(LatLon_.__name__, NN).replace('(', NN).replace(')', NN)
            self.test('quadOf', q, '77.9°S, 163.0°W, 63.1°S, 163.0°W, 63.1°S, 172.0°E, 77.9°S, 172.0°E')

        p = LatLon('66.6S', '88W')
        self.test('latlon', p.latlon, '(-66.6, -88.0)')
        self.test('philam', p.philam, '(-1.162389, -1.53589)')
#       self.test('to2ab ', p.to2ab(), p.philam)
        if LatLon is LatLon_:
            self.test('toStr', p.toStr(prec=6, kwds='test'), "66.6°S, 088.0°W, kwds='test'")
            q = p.classof(p.lat, p.lon, name='test')
            self.test('__ne__', q != p, False)
            self.test('isequalTo', q.isequalTo(p), True)  # coverage
            self.test('isequalTo', q.isequalTo(p, eps=EPS), True)

        self.test('latlonheight', p.latlonheight, '(-66.6, -88.0, 0)')
        self.test('philamheight', p.philamheight, '(-1.162389, -1.53589, 0)')

        self.test('_N_vector', p._N_vector, '(0.01386, -0.39691, -0.91775)')
        self.test('toNvector', p.toNvector().toStr(prec=5), '(0.01386, -0.39691, -0.91775)')
        self.test('toNvector', p.toNvector(Nvector=None),   '(0.01386, -0.396906, -0.917755)', known=True)

        q = p.classof(-66.6, -88)
        self.test('classof', q, p)
        try:
            t = p.others(q)
        except Exception as x:
            t = str(x)
        self.test('others', t, q)

        self.testCopy(p)


if __name__ == '__main__':  # PYCHOK internal error?

    from testRoutes import PtsFFI
    import sys

    t = Tests(__file__, __version__, points)

    p = LatLon_(0, 0)
    if hasattr(LatLon_, '__slots__'):
        try:
            t.test('LatLon_', p.__dict__, 'AttributeError')
        except AttributeError as x:
            t.test('LatLon_', x, "'LatLon_' object has no attribute '__dict__'")
        a = '__slots__'
    else:
        a = '__dict__'
    try:
        q = sys.getsizeof(p)
        t.test('sizeof(LatLon_.%s)' % (a,), q, q)
    except TypeError:  # PyPy warning/error
        pass

    pts = LatLon2psxy(PtsFFI, wrap=False)
    t.test2(pts, PtsFFI, True)

    if numpy:
        t.test('numpy.__version__', numpy.__version__, numpy.__version__)

        npy = numpy.array([(ll.lon, 0, ll.lat, 0) for ll in PtsFFI], dtype=float)
        pts = Numpy2LatLon(npy, ilat=2, ilon=0)
        t.test2(pts, npy, False)

    else:  # check abstact base class conformance
        t.test('no', 'numpy', 'numpy')

    tup = [(0, ll.lon, 0, ll.lat) for ll in PtsFFI]
    pts = Tuple2LatLon(tup, ilat=3, ilon=1)
    t.test2(pts, tup, False)
    tup = pts = None

    from pygeodesy import ellipsoidalKarney, ellipsoidalNvector, \
                          ellipsoidalVincenty, Ellipsoids, luneOf, \
                          sphericalTrigonometry, sphericalNvector
    t.test3(LatLon_)
    t.test3(ellipsoidalKarney.LatLon)
    t.test3(ellipsoidalVincenty.LatLon)
    t.test3(sphericalTrigonometry.LatLon)
    t.test3(ellipsoidalNvector.LatLon)
    t.test3(sphericalNvector.LatLon)

    if geographiclib:
        from pygeodesy.ellipsoidalKarney import areaOf as area
        q = area(luneOf(-45, 45))  # datum=WGS84
        a = Ellipsoids.WGS84.area / 4
        t.test('luneOf', q, a, fmt='%.14e', known=abs(q - a) < 0.1)  # 0.1 m^2 == 10 cm^2

    t.results()
    t.exit()
