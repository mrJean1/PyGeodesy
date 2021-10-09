
# -*- coding: utf-8 -*-

u'''Test L{iters} module.
'''

__all__ = ('Tests',)
__version__ = '21.09.30'

from base import TestsBase

from pygeodesy import NN, PointsError, PointsIter


class Tests(TestsBase):

    def testIters(self):

        Ps = PointsIter(range(8))
        for i, p in Ps.enumerate():
            pass
        self.test('i ', i, 7)
        self.test('dedup', Ps.dedup, False)

        Ps = PointsIter(range(8), loop=1)
        p0 = Ps[0]
        for i, p in Ps.enumerate(closed=True):
            pass
        self.test('i ', i, 0)
        self.test('p0', p is p0, True)
        self.test('dedup', Ps.dedup, True)  # == closed

        Ps = PointsIter(range(8))
        for p in Ps.iterate(copies=True):
            pass
        self.test('copies', Ps.copies, list(range(8)))

        r = range(8)  # list in Python 2
        Ps = PointsIter(r, loop=1)
        p0 = Ps[0]
        for i, p in Ps.enumerate(closed=True, copies=True):
            pass
        self.test('i ', i, 0)
        cs = Ps.copies
        self.test('copies', len(cs), (8 if cs is r else 9))
        self.test('p0', p is p0, True)

        t  = tuple(range(8))
        Ps = PointsIter(t)
        for p in Ps.iterate(copies=True):
            pass
        cs = Ps.copies
        self.test('copies', cs is t, True)  # original
        self.test('copies', cs, t)

        Ps = PointsIter(range(4), loop=1)
        for i, p in enumerate(Ps):
            self.test('iter', p, i + 1)
        try:
            for p in Ps:
                self.test('re-iter', p, PointsError.__name__)
        except Exception as x:
            p = repr(x).replace(',', NN)
        self.test('re-iter', p, "PointsError('points (0): too few')")


if __name__ == '__main__':

    from pygeodesy import iters

    t = Tests(__file__, __version__, iters)
    t.testIters()
    t.results()
    t.exit()
