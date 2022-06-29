
# -*- coding: utf-8 -*-

# Test some of the basics.

__all__ = ('Tests',)
__version__ = '22.06.24'

from base import TestsBase

from pygeodesy import EPS, EPS0, INF, INT0, NAN, NEG0, NINF, clips, halfs2, \
                      isclose, isfinite, isint, isint0, isneg0, isninf, isscalar, \
                      map1, property_RO, remainder, splice
from pygeodesy.basics import _xdup


class C(object):
    a = None
    b = None

    @property_RO
    def r_o(self):
        return True


class Tests(TestsBase):

    def testBasics(self):

        self.test('clips',  clips('test/testBasics.py', limit=12), 'test/t....ics.py')
        self.test('halfs2', halfs2('test/testBasics.py'), "('test/test', 'Basics.py')")

        for f, x, y, n, s in ((0,      True,  False, False, True),
                              (0.0,    True,  False, False, True),
                              (1,      True,  False, False, True),
                              (1.0,    True,  False, False, True),
                              (1e300,  True,  True,  False, True),
                              (-1e300, True,  True,  False, True),
                              (1e1234, False, False, False, True),  # == INF
                              (INF,    False, False, False, True),
                              (NAN,    False, False, False, True),
                              (NEG0,   True,  False, False, True),
                              (NINF,   False, False,  True, True)):
            t = str(f)
            self.test('isfinite(%s)'  % (t,), isfinite(f), x, nl=1)
            self.test('isint(%s)'     % (t,), isint(f, both=True), x)
            self.test('isint(%s+0.5)' % (t,), isint(f + 0.5, both=True), y)
            self.test('isninf(%s)'    % (t,), isninf(f), n)
            self.test('isscalar(%s)'  % (t,), isscalar(f), s)

        self.test('isfinite(complex)', isfinite(complex(1, 2)), True, nl=1)
        self.test('isfinite(complex)', isfinite(complex(1, NAN)), False)

        self.test('isint0(INT0)',  isint0(INT0), True, nl=1)
        self.test('isint0(False)', isint0(False), False)
        self.test('isint0(None)',  isint0(None), False)
        self.test('isint0(0)',     isint0(0), True)
        self.test('isint0(0.)',    isint0(0.), False)
        self.test('isint0(0.0)',   isint0(0, both=True), True)

        self.test('isneg0(NEG0)', isneg0(NEG0), True, nl=1)
        self.test('isneg0(0.0)',  isneg0(0.0), False)
        self.test('isneg0(INF)',  isneg0(INF), False)
        self.test('isneg0(NAN)',  isneg0(NAN), False)

        self.test('type(%s)' % ('C.r_o',), type(C.r_o).__name__, property_RO.__name__, nl=1)
        c = C()
        self.test('type(%s)' % ('c.r_o',), type(c.r_o), bool)
        self.test('c.r_o', c.r_o, True)
        try:
            c.r_o = False
            self.test('c.r_o = False', c.r_o, AttributeError.__name__)
        except AttributeError as x:
            self.test('c.r_o = False', str(x), str(x))
        except Exception as x:
            self.test('c.r_o = False', repr(x), AttributeError.__name__)

        self.test('c.a, c.b', (c.a, c.b), '(None, None)')
        d = _xdup(c, a=True, b=False)
        self.test('d.a, d.b', (d.a, d.b), '(True, False)')
        self.test('c.a, c.b', (c.a, c.b), '(None, None)')

        a, b = splice(range(10))  # PYCHOK False
        self.test('splice', (a, b), map1(type(a), (0, 2, 4, 6, 8), (1, 3, 5, 7, 9)))
        a, b, c = splice(range(10), n=3)  # PYCHOK False
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7), (2, 5, 8)))
        a, b, c = splice(range(10), n=3, fill=-1)  # PYCHOK False
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1)))
        t = tuple(splice(range(12), n=5))  # PYCHOK False
        self.test('splice', t, map1(type(t[0]), (0, 5, 10), (1, 6, 11), (2, 7), (3, 8), (4, 9)), nt=1)

        for a, b, x in (( 181,  360, '-179.0'),  # Python 3.10.4 results
                        ( 181, -360, '-179.0'),
                        ( 181,  INF,  '181.0'),
                        ( 181,  NAN,   NAN),
                        ( 181, NINF,  '181.0'),
                        (-181,  360,  '179.0'),
                        (-181, -360,  '179.0'),
                        (-181,  INF, '-181.0'),
                        (-181,  NAN,   NAN),
                        (-181, NINF, '-181.0'),
                        ( 179,  360,  '179.0'),
                        ( 179, -360,  '179.0'),
                        ( 179,  INF,  '179.0'),
                        ( 179,  NAN,   NAN),
                        ( 179, NINF,  '179.0'),
                        (-179,  360, '-179.0'),
                        (-179, -360, '-179.0'),
                        (-179,  INF, '-179.0'),
                        (-179,  NAN,   NAN),
                        (-179, NINF, '-179.0'),
                        (INF,   360,  'math domain error'),
                        (INF,  -360,  'math domain error'),
                        (INF,   INF,  'math domain error'),
                        (INF,   NAN,   NAN),
                        (INF,  NINF,  'math domain error'),
                        (NAN,   360,   NAN),
                        (NAN,  -360,   NAN),
                        (NAN,   INF,   NAN),
                        (NAN,   NAN,   NAN),
                        (NAN,  NINF,   NAN),
                        (NINF,  360,  'math domain error'),
                        (NINF, -360,  'math domain error'),
                        (NINF,  INF,  'math domain error'),
                        (NINF,  NAN,   NAN),
                        (NINF, NINF,  'math domain error')):
            try:
                r = remainder(a, b)
            except Exception as e:
                r =  str(e)
            self.test('remainder%s' % ((a, b),), r, x)

        self.test('isclose', isclose(0, EPS0), True)
        self.test('isclose', isclose(0, EPS), False)


if __name__ == '__main__':

    from pygeodesy import basics

    t = Tests(__file__, __version__, basics)
    t.testBasics()
    t.results()
    t.exit()
