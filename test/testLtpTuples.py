
# -*- coding: utf-8 -*-

u'''Test I{local tangent plane} (LTP) classes, tuples and conversions.
'''

__all__ = ('Tests',)
__version__ = '21.04.17'

from base import TestsBase

from pygeodesy import Aer, Enu, fstr, Local9Tuple, Ltp, Ned, XyzLocal
from pygeodesy.ellipsoidalKarney import Cartesian, LatLon
from pygeodesy.interns import _DOT_


def _truncate(txt):
    if len(txt) > 512:
        txt = txt[:512] + ' ...'
    return txt


# <https://www.MathWorks.com/help/map/ref/enu2geodetic.html>
Z = Ltp(46.017, 7.750, 1673, name='Zermatt')


class Tests(TestsBase):

    def testLoc(self, Loc, *args, **kwds):

        N = Loc.__name__
        z = Loc(*args, **kwds)

        t = _truncate(z.toRepr(prec=2))
        self.test(_DOT_(N, z.toRepr.__name__), t, t, nl=1)

        s = _truncate(z.toStr(prec=2))
        self.test(_DOT_(N, z.toStr.__name__), s, s)

        for C in (Aer, Enu, Ned, Xyz):  # from Loc to C and back
            n = C.__name__
            c = getattr(z.xyzLocal, 'to' + n)(C, name=n)  # Loc to C
            t = c.toStr(prec=2)
            self.test(_DOT_(N, 'xyzLocal.to' + n), t, t)

            r = getattr(c.xyzLocal, 'to' + N)(Loc, name=N)  # C back to Loc
            t = _truncate(r.toStr(prec=2))
            self.test(_DOT_(C.__name__, 'xyzLocal.to' + N), t, s)

        for C in (Cartesian, LatLon):
            n = C.__name__
            c = getattr(z.xyzLocal, 'to' + n)(C, ltp=Z, name=n)  # XXX ltp req'd
            t = _truncate(c.toStr(prec=2))
            self.test(_DOT_(N, 'to' + n), t, t)

            if Loc is Local9Tuple:
                continue

            r = c.toLocal(Loc, ltp=Z, name=N)  # XXX ltp req'd
            t = r.toStr(prec=2)
            self.test(_DOT_(C.__name__, 'toLocal ' + N), t, s)

        for a in ('azimuth', 'elevation', 'slantrange', 'groundrange',
                  'east', 'north', 'up', 'down',
                  'x', 'y', 'z', 'xyz'):  # coverage
            t = fstr(getattr(z, a), prec=3)
            self.test(_DOT_(N, a), t, t)

        if Loc is Local9Tuple:  # coverage
            for a in ('lat', 'lon', 'latlon', 'latlonheight',
                      'phi', 'lam', 'philam', 'philamheight'):
                t = fstr(getattr(z, a), prec=3)
                self.test(_DOT_(N, a), t, t)


class Xyz(XyzLocal):  # shorten name
    pass


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testLoc(Aer,  60,  40, 1000, ltp=Z)
    t.testLoc(Enu, 100, 200, 1000, ltp=Z)
    t.testLoc(Ned, 200, 100, 1000, ltp=Z)
    t.testLoc(Xyz,  10,  20, 100, ltp=Z)
    t.testLoc(Local9Tuple, 10.0, 20.0, 100.0, 46.02, 7.75, 1773.0, Z, None, None)
    t.results()
    t.exit()
