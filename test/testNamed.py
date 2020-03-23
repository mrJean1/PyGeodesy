
# -*- coding: utf-8 -*-

# Test named module.

__all__ = ('Tests',)
__version__ = '20.03.22'

from base import PyGeodesy_dir, TestsBase
from pygeodesy import Datums, named

from os import linesep

_B_    = ')}'
_C_    = '}C{'
_DICT  = 'Dict'
_LINK  = 'L{'
_TUPLE = 'Tuple'


def _mod_line(m, py):
    '''Count number of lines.
    '''
    n = py.count(linesep) + 1
    return '%s:%s' % (m, n)


class Tests(TestsBase):

    def testNamed(self, N, *args, **kwds):
        self.subtitle(named, 'ing %s%s ' % (N.__name__, args))
        n = N(*args)
        k = kwds.pop('known', False)
        self.test(n.named, n.classname, N.__name__)
        self.test(n.named, isinstance(n, N), True)
        self.test(n.named, repr(n.name), "''", known=k)
        n.name = 'Test'
        self.test(n.named, n.name, 'Test')
        self.test(n.named2, n.named2, n.classname + " 'Test'")

        for a in ('name', '_name'):  # coverage
            self.test(n.named2, getattr(n, a), 'Test')
        try:  # coverage
            self.test(n.named2, str(n), n)
        except AssertionError as x:
            self.test(n.named2, x, x)
        try:  # coverage
            self.test(n.named2, repr(n), n, f='%r', known=True)
        except AssertionError as x:
            self.test(n.named2, x, x)
        delattr(n, '_name')  # coverage
        self.test(n.named2, n.name, '')

    def testNamedDicts(self, named):
        self.subtitle(named, 'ing %s ' % ('NamedDicts',))
        self.testNamed_class(named, _DICT, '_Keys_', self._NamedDicts)

    _NamedDicts = {}  # [<name>] = 'L{<name>}C{(...)}'

    def testNamedTuples(self, named):
        self.subtitle(named, 'ing %s ' % ('NamedTuples',))
        self.testNamed_class(named, _TUPLE, '_Names_', self._NamedTuples)

    _NamedTuples = {}  # [<name>] = 'L{<name>}C{(...)}'

    def testNamed_class(self, named, _Nclass, _attr_, _Ndict):
        for n in sorted(dir(named)):
            if n.endswith(_Nclass) and n[-1 - len(_Nclass)].isdigit():
                # compare _Nattr_ and __doc__
                c = getattr(named, n)
                self.test(n, c.__name__, n)
                # check signature
                a = getattr(c, _attr_, ())
                s = 'C{(%s%s' % (', '.join(a), _B_)
                t = '%s-%s %s' % (len(a), _Nclass, s)
                d = ' '.join(c.__doc__.strip().split())
                self.test(n, t, d[:len(t)])
                # check the count
                d = n[:-len(_Nclass)]  # remove Tuple
                while not d[:1].isdigit():
                    d = d[1:]
                self.test(n, d, len(a))
                # build _Named... dict
                _Ndict[n] = 'L{%s}%s' % (n, s)

    def testNamed__doc__(self, m, py):
        self.subtitle(named, 'ing %s ' % (m,))
        for _N, _Ndict in ((_DICT,  self._NamedDicts),
                           (_TUPLE, self._NamedTuples)):
            b = max(map(len, _Ndict.keys()))
            j, _N_C_ = 0, _N + _C_
            while True:
                i = py.find(_N_C_, j)
                if i > j:
                    t.testNamed_Link(_N, i, m, py, _Ndict, b)
                    j = i + len(_N_C_)
                else:
                    break

    def testNamed_Link(self, _N, i, m, py, _Ndict, b=30):
        # check a _Named link in a __doc__ string
        L = py.rfind(_LINK, max(0, i - b), i)
        c = py.find(_B_, i)
        if 0 < L < i < c:
            m = _mod_line(m, py[:L])
            n = py[L + len(_LINK):i + len(_N)]
            t = ' '.join(py[L:c + len(_B_)].split())
            self.test(m, t, _Ndict.get(n, 'signature'))

    def testNamed_xtend(self, named):
        self.subtitle(named, 'ing %s ' % ('xtend',))
        # test extending a _NamedTuple class
        t = named.LatLon2Tuple(0, 1)
        x = t.to3Tuple(2)
        r = named.LatLon3Tuple(0, 1, 2.0)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = named.LatLon2Tuple(0, 1)
        x = t.to4Tuple(2, Datums.WGS84)
        r = named.LatLon4Tuple(0, 1, 2.0, Datums.WGS84)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = named.LatLon3Tuple(0, 1, 2)
        x = t.to4Tuple(Datums.WGS84)
        r = named.LatLon4Tuple(0, 1, 2, Datums.WGS84)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = named.PhiLam2Tuple(0, 1)
        x = t.to3Tuple(2)
        r = named.PhiLam3Tuple(0, 1, 2.0)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = named.Vector3Tuple(0, 1, 2)
        x = t.to4Tuple(4)
        r = named.Vector4Tuple(0, 1, 2, 4.0)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)


if __name__ == '__main__':

    from glob import glob
    import os.path as os_path

    from pygeodesy import ecef, elliptic, frechet, hausdorff, \
                          Transforms

    t = Tests(__file__, __version__)
    t.testNamed(named._Named)
    t.testNamed(named._NamedBase)
    t.testNamed(named._NamedDict)
    t.testNamed(named._NamedEnum, 'Test', known=True)
    t.testNamed(named._NamedEnumItem)
    t.testNamed(named._NamedInt, 0)
    t.testNamed(named._NamedStr, '')
    t.testNamed(named.LatLon2Tuple, 0, 0)  # _NamedTuple

    # find all _NamedDict and _NamedTuple (sub)classes
    t.testNamedDicts(named)
    t.testNamedTuples(named)
    t.testNamedTuples(ecef)
    t.testNamedTuples(elliptic)
    t.testNamedTuples(frechet)
    t.testNamedTuples(hausdorff)

    # test __doc__ strings in all pygeodesy modules
    for m in glob(os_path.join(PyGeodesy_dir, 'pygeodesy', '*.py')):
        with open(m, 'rb') as f:
            py = f.read()
            if isinstance(py, bytes):  # Python 3+
                py = py.decode('utf-8')
            t.testNamed__doc__(os_path.basename(m), py)

    t.testNamed_xtend(named)

    t.subtitle(named, 'ing %s ' % ('unregister',))
    for n, d in tuple(Datums.items()):
        Datums.unregister(d)  # coverage _NamedEnum.unregister
        t.test('unregister(%s)' % (n,), getattr(Datums, n, None), None)
    t.test('Datums', len(Datums), 0)

    for n, e in tuple(Transforms.items()):
        e.unregister()  # coverage _NamedEnumItem.unregister
        t.test(n + '.unregister', getattr(Transforms, n, None), None)
    t.test('Transforms', len(Transforms), 0)

    t.subtitle(named, 'ing %s ' % ('coverage',))
    Nd = named.Neighbors8Dict  # coverage
    nd = Nd(**dict((t, t) for t in (Nd._Keys_ + ('name',))))
    t.test('nd.name', nd.name, 'name')
    t.test('nd.named', nd.named, 'name')
    del nd.name
    t.test('nd.named', nd.named, Nd.__name__)

    nd.name = 'test'
    t.test('nd.name', nd.name, 'test')

    nd.test = 'test'
    t.test('nd.test', nd.test, 'test')
    del nd.test
    t.test('nd.test', getattr(nd, 'test', None), None)

    t.test('nd.classnaming', nd.classnaming, False)
    t.test('nd.classname', nd.classname, nd.classname)
    t.test('nd.named2', nd.named2, nd.named2)

    nd.classnaming = True
    t.test('nd.classnaming', nd.classnaming, True)
    t.test('nd.classname', nd.classname, nd.classname)
    t.test('nd.named2', nd.named2, nd.named2)

    t.test('classnaming', named.classnaming(True), False)
    t.test('classnaming', named.classnaming(False), True)

    t.results()
    t.exit()
