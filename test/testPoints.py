
# -*- coding: utf-8 -*-

# Test the simplify functions.

__all__ = ('Tests',)
__version__ = '17.09.14'

from base import TestsBase

from pygeodesy import EPS, LatLon2psxy, Numpy2LatLon, classname, points

try:
    from collections import Sequence
except ImportError:
    Sequence = None

try:
    import numpy
except ImportError:
    numpy = None


class Tests(TestsBase):

    def test2(self, pts, npy, psxy=False):

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
        _test('len',    len(pts), len(npy))
        _test('iter',   len(tuple(iter(pts))), len(npy))
        if hasattr(npy, 'shape'):
            _test('shape',  npy.shape, pts.shape)
        _test('slice1', len(pts[:n]), n)
        _test('slice2', type(pts[1:n:2]), type(pts))
        _test('slice3', pts[1:n][0], pts[1])
        _test('strepr', str(pts), repr(pts))
        if hasattr(pts, 'subset'):
            _test('subset', type(pts.subset(range(n))), type(npy))  # , nt=1)

        for i, p in ((10, (52.224006, -0.707747)),
                     (20, (52.232688, -0.714608)),
                     (30, (52.234375, -0.714348)),
                     (40, (52.237239, -0.712557)),
                     (50, (52.24023,  -0.709919)),
                     (60, (52.240745, -0.707042))):

            if psxy:  # flip to x, y tuple
                p = tuple(reversed(p))

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

        i = len(pts)
        for p in reversed(pts):
            i -= 1
            _test('reversed[%s]' % (i,), p, pts[i])


if __name__ == '__main__':  # PYCHOK internal error?

    from testRoutes import PtsFFI

    t = Tests(__file__, __version__, points)

    pts = LatLon2psxy(PtsFFI, wrap=False)
    t.test2(pts, PtsFFI, True)

    if numpy:
        t.test('numpy.__version__', numpy.__version__, numpy.__version__)

        npy = numpy.array([(ll.lon, 0, ll.lat, 0) for ll in PtsFFI], dtype=float)
        pts = Numpy2LatLon(npy, ilat=2, ilon=0)

        t.test2(pts, npy, False)

    else:  # check abstact base class conformance
        t.test('no module', 'numpy', 'numpy')

    t.results(nl=0)
    t.exit()
