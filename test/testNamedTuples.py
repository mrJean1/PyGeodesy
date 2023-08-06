
# -*- coding: utf-8 -*-

# Test L{namedTuples} module.

__all__ = ('Tests',)
__version__ = '23.08.05'

from bases import TestsBase
from pygeodesy import FIx, issubclassof
from pygeodesy.albers import _Ks
from pygeodesy.frechet import Frechet6Tuple
from pygeodesy.fsums import _Float_Int
from pygeodesy.hausdorff import Hausdorff6Tuple
from pygeodesy.interns import _DOT_
from pygeodesy.karney import _GTuple
from pygeodesy.named import _NamedTuple, _Pass
# from pygeodesy.rhumbBase import _Lat as rLat, _Lon as rLon
from pygeodesy.triaxials import _NamedTupleTo
from pygeodesy.units import _NamedUnit

_Units_ = '_Units_'
_Xcepts = (FIx, _Float_Int, _Ks, _Pass)  # rLat, rLon


class Tests(TestsBase):

    def testNamedTuple(self, T, *args):
        m = T.__module__
        n = _DOT_(T.__name__, _Units_)
        t = T(*args)

        U = getattr(T, _Units_, ())
        if len(U) != len(t):  # check _Units_ len
            e = '%s = %r' % (n, tuple(u.__name__ for u in U))
            self.test(m, e, len(t))

        for i, u in enumerate(U):  # check _Units_ types
            if not (callable(u) and (u in _Xcepts or
                                     issubclassof(u, _NamedUnit))):
                e = '%s[%s] %r' % (n, i, u)
                self.test(m, e, callable.__name__)

        u = t.toUnits()  # check _Units_ sample
        ru = repr(u)
        rt = repr(t)
        x = ru if rt == ru.replace("'0", '0').replace("5'", '5') \
                          .replace('True', '0.5').replace('0,', '0.5,') \
                          .replace('0)', '0.5)') else rt
        self.test(m, ru, x)

        u = ('%s=%s' % (n, u.__name__) for (n, u) in zip(T._Names_, U))
        u = '%s(%s)' % (T.__name__, ', '.join(u))
        self.test(m, u, u)  # the items as name=units

        c = _DOT_(m, t.__class__.__name__)
        for n in T._Names_:  # coverage
            x = str(getattr(t, n))
            self.test(_DOT_(c, n), x, x)

        self.testValidated(T, True)
        self.testValidated(t.__class__, True)

    def testNamedTuples(self):
        for T in self.pygeodesy_classes(Base=_NamedTuple):
            if T not in (_GTuple, _NamedTuple, _NamedTupleTo):
                t = (0.5,) * len(T._Names_)  # some sample value
                if T in (Frechet6Tuple, Hausdorff6Tuple):
                    t = t[1:] + ('test',)
                self.testNamedTuple(T, *t)

        self.testValidated(_NamedTuple, False)

    def testValidated(self, T, x):
        n = _DOT_(T.__module__, _DOT_(T.__name__, '_validated'))
        self.test(n , T._validated, x)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testNamedTuples()
    t.results()
    t.exit()
