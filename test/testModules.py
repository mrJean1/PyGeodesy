
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '17.06.23'

from base import TestsBase, type2str


class Tests(TestsBase):

    def testModule(self, m, name=''):
        # check that __all__ names exist in module m
        self.subtitle(m, 'Module')

        n_ = (name or m.__name__).split('.')[-1] + '.'
        for a in sorted(m.__all__):
            n = n_ + a + type2str(m, a)
            o = getattr(m, a, None)
            t = getattr(o, '__module__', None)
            if t and t != m.__name__:
                n = '%s (%s)' % (n, t)
            self.test(n, hasattr(m, a), True)


if __name__ == '__main__':

    from pygeodesy import datum, dms, \
                          ellipsoidalNvector, ellipsoidalVincenty, \
                          lcc, mgrs, nvector, osgr, simplify, \
                          sphericalNvector, sphericalTrigonometry, \
                          vector3d, utm, utils  # PYCHOK expected
    import pygeodesy  # PYCHOK expected

    t = Tests(__file__, __version__)

    t.testModule(pygeodesy, 'pygeodesy')
    for m in (datum, dms,
              ellipsoidalNvector, ellipsoidalVincenty,
              lcc, mgrs, nvector, osgr, simplify,
              sphericalNvector, sphericalTrigonometry,
              vector3d, utm, utils):
        t.testModule(m)

    t.results()
    t.exit()
