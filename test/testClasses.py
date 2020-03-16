
# -*- coding: utf-8 -*-

# Test class attributes and mro's, especially since
# the class hierarchy is non-trivial for certain
# classes, like LatLon.

__all__ = ('Tests',)
__version__ = '20.03.14'

from inspect import isclass
from os.path import basename

from base import TestsBase, type2str


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
        for a, m in sorted(attrs.items()):
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
            if isclass(C) and not (a.endswith('Error')
                                or a in ('property_RO', 'SciPyWarning')):
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

    from pygeodesy import ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty, \
                          nvector, \
                          sphericalNvector, sphericalTrigonometry, \
                          vector3d  # PYCHOK expected

    t = Tests(__file__, __version__)

    # check Cartesian attributes and mro's
    t.testCartesianAttrs(ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty,
                         sphericalNvector, sphericalTrigonometry)

    # check LatLon attributes and mro's
    t.testLatLonAttrs(ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty,
                      sphericalNvector, sphericalTrigonometry)

    # check vector attributes and mro's
    t.testVectorAttrs(ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty,
                      sphericalNvector, sphericalTrigonometry,
                      nvector, vector3d)

    import pygeodesy  # PYCHOK re-import
    t.testCopyAttr(pygeodesy)

    t.results()
    t.exit()
