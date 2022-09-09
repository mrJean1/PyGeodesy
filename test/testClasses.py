
# -*- coding: utf-8 -*-

# Test class attributes and mro's, especially since
# the class hierarchy is non-trivial for certain
# classes, like LatLon.

__all__ = ('Tests',)
__version__ = '22.09.02'

from base import GeodSolve, TestsBase, type2str

from pygeodesy import itemsorted, Property, Property_RO, property_RO, \
                      SciPyWarning, Str_

from inspect import isclass
from os.path import basename

_No_Copy_OK = set((Property, Property_RO, property_RO, SciPyWarning, Str_))


class Tests(TestsBase):

    def _subtitle(self, classname, testing, unused):
        # can't overload method TestsBase.subtitle
        self.printf('test%s%s(%s)', classname, testing, __version__, nl=1)

    def testAttrs(self, classname, modules, *args, **kwds):
        self._subtitle(classname, 'Attrs', modules)
        attrs = {}
        for m in modules:
            C = getattr(m, classname, None)
            if C:
                i = C(*args, **kwds)
                n = basename(m.__name__)
                for a in dir(i):
                    if not a.startswith('_'):
                        a += type2str(C, a)
                        attrs[a] = attrs.get(a, ()) + (n,)
        for a, m in itemsorted(attrs):
            m = ', '.join(sorted(m))
            self.test(a, m, m)  # passes always

    def testMro(self, classname, modules):
        self._subtitle(classname, 'Mro', modules)
        for m in modules:
            C = getattr(m, classname, None)
            if C:
                c = ', '.join(str(c)[8:-2] for c in C.mro()[:-1])
                self.test(m.__name__, c, c)  # passes always

    def testCartesianAttrs(self, *modules):
        self.testAttrs('Cartesian', modules, 0, 0, 0)
        self.testMro(  'Cartesian', modules)

    def testCopyAttr(self, package):
        self._subtitle('Copy', 'Attr', package)
        for a in sorted(package.__all__):
            C = getattr(package, a)
            if isclass(C) and not (C in _No_Copy_OK or a.endswith('Error')):
                c = getattr(C, 'copy', None) or getattr(C, 'fcopy', None)
                t = 'copy' if callable(c) else 'missing'
                self.test(C.__name__, t, 'copy')

    def testLatLonAttrs(self, *modules):
        self.testAttrs('LatLon', modules, 0, 0)
        self.testMro(  'LatLon', modules)

    def testVectorAttrs(self, *modules):
        self.testAttrs('Nvector', modules, 0, 0, 0, h=0)
        self.testMro(  'Nvector', modules)
        self.testAttrs('Vector3d', modules, 0, 0, 0)
        self.testMro(  'Vector3d', modules)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalKarney, \
                          ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry, \
                          vector3d, nvectorBase  # DEPRECATED nvector

    t = Tests(__file__, __version__)

    ms = (sphericalNvector, sphericalTrigonometry,
          ellipsoidalNvector, ellipsoidalVincenty,
          ellipsoidalKarney, ellipsoidalExact)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        ms += (ellipsoidalGeodSolve,)

    # check Cartesian attributes and mro's
    t.testCartesianAttrs(*ms)

    # check LatLon attributes and mro's
    t.testLatLonAttrs(*ms)

    # check vector attributes and mro's
    t.testVectorAttrs(nvectorBase, vector3d, *ms)  # DEPRECATED nvector

    import pygeodesy  # PYCHOK re-import
    t.testCopyAttr(pygeodesy)

    t.results()
    t.exit()
