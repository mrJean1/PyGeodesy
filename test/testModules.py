
# -*- coding: utf-8 -*-

# Test modules and attributes.

__all__ = ('Tests',)
__version__ = '23.09.13'

from bases import isPyPy, TestsBase, type2str
from pygeodesy.interns import _DOT_


class Tests(TestsBase):

    def testModule(self, m, name=''):
        # check that __all__ names exist in module m
        self.subtitle(m, 'Module')

        m_ = (name or m.__name__).split(_DOT_)[-1] + _DOT_
        for a in sorted(m.__all__):
            n = m_ + a + type2str(m, a)
            o = getattr(m, a, None)
            t = getattr(o, '__module__', None)
            if t and t != m.__name__:
                n = '%s (%s)' % (n, t)
            self.test(n, hasattr(m, a), True)


if __name__ == '__main__':

    import pygeodesy

    t = Tests(__file__, __version__)
    if pygeodesy.__all__:  # and pygeodesy._init__all__
        pys = set(pygeodesy.__all__)

        from inspect import ismodule

        t.testModule(pygeodesy, 'pygeodesy')
        for a in sorted(pys):
            m = getattr(pygeodesy, a)
            if ismodule(m):
                t.testModule(m)

        if not isPyPy:  # XXX because next line throws ...
            # AttributeError: 'method' object has no attribute '__module__'
            # at least when testing with PyPy on TravisCI

            # keyword.iskeyword.__module__ == None
            pygeodesy.iskeyword.__module__ = 'keyword'

            from pygeodesy.interns import _sub_packages

            # check module for public functions, etc.
            t.subtitle(pygeodesy, 'Public')
            for a in sorted(pys):
                f = getattr(pygeodesy, a)
                m = getattr(f, '__module__', _DOT_).split(_DOT_)[:2][-1]
                if m and m not in ('math', 'pygeodesy'):
                    n = a + type2str(pygeodesy, a)
                    t.test(n, m in  pys or
                              m in _sub_packages or
                              a in ('freduce' , 'isclass', 'iskeyword'), True)

    else:
        t.skip('pygeodesy', 2805)
    t.results()
    t.exit()
