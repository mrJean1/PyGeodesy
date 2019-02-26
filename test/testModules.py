
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '19.02.21'

from base import TestsBase, type2str


class Tests(TestsBase):

    def testModule(self, m, name=''):
        # check that __all__ names exist in module m
        self.subtitle(m, 'Module')

        m_ = (name or m.__name__).split('.')[-1] + '.'
        for a in sorted(m.__all__):
            n = m_ + a + type2str(m, a)
            o = getattr(m, a, None)
            t = getattr(o, '__module__', None)
            if t and t != m.__name__:
                n = '%s (%s)' % (n, t)
            self.test(n, hasattr(m, a), True)


if __name__ == '__main__':

    import pygeodesy  # PYCHOK expected
    from inspect import ismodule

    t = Tests(__file__, __version__)

    t.testModule(pygeodesy, 'pygeodesy')
    for a in sorted(pygeodesy.__all__):
        m = getattr(pygeodesy, a)
        if ismodule(m):
            t.testModule(m)

    # check module for public functions, etc.
    t.subtitle(pygeodesy, 'Public')
    for a in sorted(pygeodesy.__all__):
        f = getattr(pygeodesy, a)
        m = getattr(f, '__module__', '.').split('.')[-1]
        if m and m not in ('bases', 'math', 'pygeodesy'):
            n = a + type2str(pygeodesy, a)
            t.test(n, m in pygeodesy.__all__ or a == 'freduce', True)

    t.results()
    t.exit()
