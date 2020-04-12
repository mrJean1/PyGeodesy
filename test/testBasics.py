
# -*- coding: utf-8 -*-

# Test some of the basics.

__all__ = ('Tests',)
__version__ = '20.04.11'

from base import TestsBase

from pygeodesy import INF, NAN, NEG0, clips, halfs2, \
                      isint, isfinite, isneg0, isscalar, \
                      LenError, LimitError, limiterrors, map1, \
                      property_RO, RangeError, rangerrors, splice


class C(object):
    @property_RO
    def r_o(self):
        return True


class Tests(TestsBase):

    def testBasics(self):

        self.test('clips',  clips('test/testBasics.py', limit=12), 'test/t....ics.py')
        self.test('halfs2', halfs2('test/testBasics.py'), "('test/test', 'Basics.py')")

        for f, x, y, s in ((0,      True,  False, True),
                           (0.0,    True,  False, True),
                           (1,      True,  False, True),
                           (1.0,    True,  False, True),
                           (1e300,  True,  True,  True),
                           (-1e300, True,  True,  True),
                           (1e1234, False, False, True),  # == INF
                           (INF,    False, False, True),
                           (NAN,    False, False, True),
                           (NEG0,   True,  False, True)):
            t = str(f)
            self.test('isfinite(%s)'  % (t,), isfinite(f), x)
            self.test('isint(%s)'     % (t,), isint(f, both=True), x)
            self.test('isint(%s+0.5)' % (t,), isint(f + 0.5, both=True), y)
            self.test('isscalar(%s)'  % (t,), isscalar(f), s)

        self.test('isneg0(NEG0)', isneg0(NEG0), True)
        self.test('isneg0(0.0)',  isneg0(0.0), False)
        self.test('isneg0(INF)',  isneg0(INF), False)
        self.test('isneg0(NAN)',  isneg0(NAN), False)

        self.test('limiterrors', limiterrors(False), True)
        self.test('limiterrors', limiterrors(True), False)
        self.test('rangerrors',  rangerrors(False),  True)
        self.test('rangerrors',  rangerrors(True),  False)

        self.test('type(%s)' % ('C.r_o',), type(C.r_o).__name__, property_RO.__name__)
        c = C()
        self.test('type(%s)' % ('c.r_o',), type(c.r_o), bool)
        self.test('c.r_o', c.r_o, True)
        try:
            c.r_o = False
            self.test('c.r_o = False', c.r_o, AttributeError)
        except AttributeError as x:
            x = repr(x).split()
            x = ' '.join(x[:2] + ['....'] + x[-2:])
            x = x.replace(',', '')
            self.test('c.r_o = False', x, "AttributeError('immutable property: .... = False')")
        except Exception as x:
            self.test('c.r_o = False', repr(x), AttributeError)

        a, b = splice(range(10))  # PYCHOK False
        self.test('splice', (a, b), map1(type(a), (0, 2, 4, 6, 8), (1, 3, 5, 7, 9)))
        a, b, c = splice(range(10), n=3)  # PYCHOK False
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7), (2, 5, 8)))
        a, b, c = splice(range(10), n=3, fill=-1)  # PYCHOK False
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1)))
        t = tuple(splice(range(12), n=5))  # PYCHOK False
        self.test('splice', t, map1(type(t[0]), (0, 5, 10), (1, 6, 11), (2, 7), (3, 8), (4, 9)))

    def testErrors(self, InvalidError, IsnotError):

        e = InvalidError(zero=1)
        self.test(InvalidError.__name__, e, 'zero invalid: 1')
        self.test(InvalidError.__name__, repr(e).replace(',', ''), "ValueError('zero invalid: 1')")
        e = InvalidError(zero=2, Error=RangeError)
        self.test(InvalidError.__name__, e, 'zero invalid: 2')
        self.test(InvalidError.__name__, repr(e).replace(',', ''), "RangeError('zero invalid: 2')")

        e = IsnotError(int.__name__, float.__name__, _None=None)
        self.test(IsnotError.__name__, e, '_None is not an int or float: None')
        self.test(IsnotError.__name__, repr(e).replace(',', ''), "TypeError('_None is not an int or float: None')")

        e = IsnotError('scalar', _None=None)
        self.test(IsnotError.__name__, e, '_None is not scalar: None')
        self.test(IsnotError.__name__, repr(e).replace(',', ''), "TypeError('_None is not scalar: None')")
        e = IsnotError('scalar', Error=LimitError, _None=None)
        self.test(IsnotError.__name__, e, '_None is not scalar: None')
        self.test(IsnotError.__name__, repr(e).replace(',', ''), "LimitError('_None is not scalar: None')")

        try:
            raise LenError(LenError, a=1, b=2, c=3, d=4)
            self.test(LenError.__name__, None, LenError)
        except ValueError as x:
            self.test(LenError.__name__, str(x), 'LenError(a, b, c, d) len: 1 vs 2 vs 3 vs 4')

    def testFloat(self, _float):
        n = _float.__name__

        t.test(n, _float(12),   '12.0')
        t.test(n, _float(f=12), '12.0')
        try:
            t.test(n, _float('X'), ValueError)
        except ValueError as x:
            t.test(n, str(x), "invalid, single float('X')")

        try:
            t.test(n, _float(12, f=12), ValueError)
        except ValueError as x:
            t.test(n, str(x), 'invalid, single float(12, f=12)')

        try:
            t.test(n, _float('X', Error=TypeError), TypeError)
        except TypeError as x:
            t.test(n, x.__class__, TypeError)


if __name__ == '__main__':

    from pygeodesy import basics

    t = Tests(__file__, __version__, basics)
    t.testBasics()
    t.testErrors(basics.InvalidError, basics.IsnotError)
    t.testFloat(basics._float)
    t.results()
    t.exit()
