
# -*- coding: utf-8 -*-

# Test named module.

__all__ = ('Tests',)
__version__ = '21.12.20'

from base import TestsBase
from pygeodesy import geohash, Datum, Datums, named, \
                      namedTuples, nameof, ub2str

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

    def testBases(self):
        self.subtitle(named, 'ing %s ' % ('bases',))

        nd = named._NamedDict({'1': 1, '2': 2}, name='test')
        self.test('nd.dict', nd.toRepr(), 'test(1=1, 2=2)')
        self.test('nd.name', nd.name, 'test')
        nd = named._NamedDict({'1': 1, '2': 2, 'name': 'test'})
        self.test('nd.dict', nd.toRepr(), 'test(1=1, 2=2)')
        self.test('nd.name', nd.name, 'test')
        nd = named._NamedDict(one=1, two=2, name='test')
        self.test('nd.kwds', nd.toRepr(), 'test(one=1, two=2)')
        self.test('nd.name', nd.name, 'test')
        nd = named._NamedDict({'1': 1, '2': 2, 'name': 'test'}, name='kwds')
        self.test('nd.dict', nd.toRepr(), 'test(1=1, 2=2)')
        self.test('nd.name', nd.name, 'test')
        nd = named._NamedDict([('1', 1), ('2', 2), ('name', 'test')], name='kwds')
        self.test('nd.list', nd.toRepr(), 'test(1=1, 2=2)')
        self.test('nd.name', nd.name, 'test')
        nd.update(dict(name='kwds'))
        self.test('nd.updated', nd.toRepr(), "test(1=1, 2=2, name='kwds')")
        self.test('nd.name', nd.name, 'test')

    def testNamed(self, N, *args, **kwds):
        self.subtitle(named, 'ing %s%s ' % (N.__name__, args))
        k = kwds.pop('known', False)
        n = N(*args)
        self.test(n.named, n.classname, N.__name__)
        self.test(n.named, isinstance(n, N), True)
        self.test(n.named, repr(n.name), "''", known=k)
        n.rename('Test')
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

    def testNamedDicts(self):
        self.subtitle(named, 'ing %s ' % ('NamedDicts',))
        self.testNamed_classes(named._NamedDict, _DICT, '_Keys_', self._NamedDicts)

    _NamedDicts = {}  # [<name>] = 'L{<name>}C{(...)}'

    def testNamedTuples(self):
        self.subtitle(namedTuples, 'ing %s ' % ('NamedTuples',))
        self.testNamed_classes(named._NamedTuple, _TUPLE, '_Names_', self._NamedTuples)

    _NamedTuples = {}  # [<name>] = 'L{<name>}C{(...)}'

    def testNamed_classes(self, _Nbase, _Nclass, _attr_, _Ndict):
        for c in self.pygeodesy_classes(Base=_Nbase):
            n = c.__name__
            if c is not _Nbase and n[-1 - len(_Nclass)].isdigit():
                # compare _Nattr_ and __doc__
                self.test(n, c.__name__, n)
                # check signature
                a = getattr(c, _attr_, ())
                s = 'C{(%s%s' % (', '.join(a), _B_)
                t = '%s-%s %s' % (len(a), _Nclass, s)
                d = ' '.join(c.__doc__.strip().split())
                self.test(n, d[:len(t)].strip(), t)
                # check the count
                d = n[:-len(_Nclass)]
                while d and not d[:1].isdigit():
                    d = d[1:]
                if d:
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

    def testNamed_xtend(self, namedTuples):
        self.subtitle(namedTuples, 'ing %s ' % ('xtend',))
        # test extending a _NamedTuple class
        t = namedTuples.LatLon2Tuple(0, 1)
        x = t.to3Tuple(2)
        r = namedTuples.LatLon3Tuple(0, 1, 2)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = namedTuples.LatLon2Tuple(0, 1)
        x = t.to4Tuple(2, Datums.WGS84)
        r = namedTuples.LatLon4Tuple(0, 1, 2, Datums.WGS84)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = namedTuples.LatLon3Tuple(0, 1, 2)
        x = t.to4Tuple(Datums.WGS84)
        r = namedTuples.LatLon4Tuple(0, 1, 2, Datums.WGS84)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = namedTuples.PhiLam2Tuple(0, 1)
        x = t.to3Tuple(2)
        r = namedTuples.PhiLam3Tuple(0, 1, 2)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

        t = namedTuples.Vector3Tuple(0, 1, 2)
        x = t.to4Tuple(4)
        r = namedTuples.Vector4Tuple(0, 1, 2, 4)
        self.test(repr(t), x, r)
        self.test(repr(t), x.__class__, r.__class__)

    def testUnregister(self):
        self.subtitle(named, 'ing %s ' % ('unregister',))
        from pygeodesy import Conics, Ellipsoid, Ellipsoids, RefFrames, Transforms

        self.test('Conics', len(Conics), 1, known=True)
        t = tuple(Conics.items(all=True))
        self.test('Conics', len(Conics), 8)
        for n, c in t:
            c.unregister()  # coverage _NamedEnum.unregister
            self.test('Conics.' + n + '.unregister', getattr(Conics, n, None), None)
        self.test('Conics', len(Conics), 0)

        self.test('Datums', len(Datums), 6, known=True)
        t = tuple(Datums.items(all=True))
        self.test('Datums', len(Datums), 18)
        for n, d in t:
            Datums.unregister(d)  # coverage _NamedEnum.unregister
            self.test('Datums.unregister(%s)' % (n,), getattr(Datums, n, None), None)
        self.test('Datums', len(Datums), 0)

        self.test('RefFrames', len(RefFrames), 0, known=True)
        t = tuple(RefFrames.items(all=True))
        self.test('RefFrames', len(RefFrames), 15)
        for n, r in t:
            r.unregister()  # coverage _NamedEnum.unregister
            self.test('RefFrames.' + n + '.unregister', getattr(RefFrames, n, None), None)
        self.test('RefFrames', len(RefFrames), 0)

        self.test('Ellipsoids', len(Ellipsoids), 12, known=True)
        t = tuple(Ellipsoids.items(all=True))
        self.test('Ellipsoids', len(Ellipsoids), 41)
        for n, e in t:
            e.unregister()  # coverage _NamedEnum.unregister
            self.test('Ellipsoids.' + n + '.unregister', getattr(Ellipsoids, n, None), None)
        self.test('Ellipsoids', len(Ellipsoids), 0)

        self.test('Transforms', len(Transforms), 15, known=True)
        t = tuple(Transforms.items(all=True))
        self.test('Transforms', len(Transforms), 18)
        for n, x in t:
            x.unregister()  # coverage _NamedEnumItem.unregister
            self.test('Transforms.' + n + '.unregister', getattr(Transforms, n, None), None)
        self.test('Transforms', len(Transforms), 0)

        self.test('nameof', nameof(Ellipsoid.KsOrder), 'KsOrder')


if __name__ == '__main__':

    from base import isiOS

    t = Tests(__file__, __version__)

    t.testNamed(named._Named)
    t.testNamed(named._NamedBase)
    t.testNamed(named._NamedDict)
    t.testNamed(named._NamedEnum, Datum, known=True)
    t.testNamed(named._NamedEnumItem)
    t.testNamed(namedTuples.LatLon2Tuple, 0, 0)  # _NamedTuple

    # find _NamedDict and _NamedTuple (sub)classes
    # defined in all pygeodesy modules
    t.testNamedDicts()
    t.testNamedTuples()

    # test __doc__ strings in all pygeodesy modules
    for n, m in t.pygeodesy_names2():
        with open(n, 'rb') as f:
            py = ub2str(f.read())
            t.testNamed__doc__(m, py)

    t.testNamed_xtend(namedTuples)

    t.testBases()

    t.subtitle(named, 'ing %s ' % ('coverage',))
    Nd = geohash.Neighbors8Dict  # coverage
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

    if not isiOS:  # Pythonista runs in single process
        t.testUnregister()

    t.results()
    t.exit()
