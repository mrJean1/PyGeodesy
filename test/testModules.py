
# -*- coding: utf-8 -*-

# Test modules and attributes.

__all__ = ('Tests',)
__version__ = '23.11.26'

from bases import isPyPy, TestsBase, type2str
from pygeodesy.interns import _DOT_, _headof, _tailof


class Tests(TestsBase):

    def testModule(self, m, name=''):
        # check that __all__ names exist in module m
        self.subtitle(m, 'Module')

        m_ = _tailof(name or m.__name__) + _DOT_
        for a in sorted(m.__all__):
            n = m_ + a + type2str(m, a)
            o =  getattr(m, a, None)
            t = _parent(o, None)
            if t and t != m.__name__:
                n = '%s (%s)' % (n, t)
            self.test(n, hasattr(m, a), True)


def _parent(obj, dflt):
    try:  # module and package
        return getattr(obj, '__module__', dflt)
    except (AttributeError, ImportError):
        return getattr(obj, '__package__', dflt)


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
                f =  getattr(pygeodesy, a)
                m = _headof(_headof(_parent(f, _DOT_)))
                if m and m not in ('math', 'pygeodesy'):
                    n = a + type2str(pygeodesy, a)
                    t.test(n, m in  pys or
                              m in _sub_packages or
                              a in ('freduce' , 'isclass', 'iskeyword'), True)

    else:
        t.skip('pygeodesy', 2805)
    t.results()
    t.exit()
