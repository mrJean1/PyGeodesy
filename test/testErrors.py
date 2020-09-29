
# -*- coding: utf-8 -*-

# Test some of the errors.

__all__ = ('Tests',)
__version__ = '20.09.27'

from base import isPython3, TestsBase

from pygeodesy import crosserrors, exception_chaining, LenError, \
                      LimitError, limiterrors, \
                      RangeError, rangerrors

from os import environ


class Tests(TestsBase):

    def testError(self, Error):
        try:
            raise Error('name', 'value', txt='test1 txt')
        except Error as x:
            t = str(x)
            self.test(Error.__name__, t, t)
        try:
            raise Error(txt='test2 txt')  # coverage
        except Error as x:
            t = str(x)
            self.test(Error.__name__, t, t)

    def testErrors(self, InvalidError, IsnotError):

        e = InvalidError(zero=1)
        self.test(InvalidError.__name__, e, 'zero (1): invalid')
        self.test(InvalidError.__name__, repr(e).replace(',)', ')'), "ValueError('zero (1): invalid')")
        e = InvalidError(zero=1, one=2, txt='outside')
        self.test(InvalidError.__name__, e, 'one (2) or zero (1): outside')
        self.test(InvalidError.__name__, repr(e).replace(',)', ')'), "ValueError('one (2) or zero (1): outside')")
        e = InvalidError(zero=1, Error=RangeError, one=2, txt='outside')
        self.test(InvalidError.__name__, e, 'one (2) or zero (1): outside')
        self.test(InvalidError.__name__, repr(e).replace(',)', ')'), "RangeError('one (2) or zero (1): outside')")
        e = InvalidError(two=2, Error=ValueError)  # coverage

        e = IsnotError(int.__name__, float.__name__, _None=None)
        self.test(IsnotError.__name__, e, '_None (None) not an int or float')
        self.test(IsnotError.__name__, repr(e).replace(',)', ')'), "TypeError('_None (None) not an int or float')")

        e = IsnotError('scalar', _None=None)
        self.test(IsnotError.__name__, e, '_None (None) not scalar')
        self.test(IsnotError.__name__, repr(e).replace(',)', ')'), "TypeError('_None (None) not scalar')")
        e = IsnotError('scalar', Error=LimitError, _None=None)
        self.test(IsnotError.__name__, e, '_None (None) not scalar: invalid')
        self.test(IsnotError.__name__, repr(e).replace(',)', ')'), "LimitError('_None (None) not scalar: invalid')")

        try:
            raise LenError(LenError, a=1, b=2, c=3, d=4)
            self.test(LenError.__name__, None, LenError.__name__)
        except ValueError as x:
            self.test(LenError.__name__, str(x), 'LenError(a, b, c, d) len 1 vs 2 vs 3 vs 4: invalid')

        self.test(crosserrors.__name__, crosserrors(False), True)
        self.test(crosserrors.__name__, crosserrors(True), False)
        self.test(limiterrors.__name__, limiterrors(False), True)
        self.test(limiterrors.__name__, limiterrors(True), False)
        self.test(rangerrors.__name__,  rangerrors(False),  True)
        self.test(rangerrors.__name__,  rangerrors(True),  False)

        x = True if environ.get('PYGEODESY_EXCEPTION_CHAINING', None) else False
        self.test(exception_chaining.__name__, exception_chaining(), x if isPython3 else None)
        self.test(exception_chaining.__name__, exception_chaining(RangeError()), None)
        self.test(exception_chaining.__name__, exception_chaining(TypeError()), None)

    def testKwds(self, xkwds):
        self.test(xkwds.__name__, xkwds({}, test='test1'), 'test1')
        self.test(xkwds.__name__, xkwds({'test': 'test2'}, test='test3'), 'test2')
        try:
            x = AssertionError.__name__
            t = xkwds({})
        except AssertionError as a:
            t = x = str(a)
        self.test(xkwds.__name__, t, x)
        try:
            x = AssertionError.__name__
            t = xkwds({}, n1='d1', n2='d2')
        except AssertionError as a:
            t = x = str(a)
        self.test(xkwds.__name__, t, x)


if __name__ == '__main__':

    from pygeodesy import errors, ClipError, CrossError, CSSError, EcefError, \
                                  EllipticError, EPSGError, ETMError, FrechetError, \
                                  GARSError, GeohashError, GeoidError, HausdorffError, \
                                  HeightError, LazyImportError, LCCError, MGRSError, \
                                  OSGRError, PGMError, PointsError, RefFrameError, \
                                  SciPyError, SciPyWarning, TRFError, UnitError, \
                                  UPSError, UTMError, UTMUPSError, VectorError, \
                                  VincentyError, WebMercatorError, WGRSError  # ClipError

    t = Tests(__file__, __version__, errors)
    for E in (errors._AssertionError, errors._AttributeError, errors._IndexError,
              errors.LimitError,      errors._NameError,      errors.ParseError,
              errors._TypeError,
              ClipError, CrossError, CSSError, EcefError, EllipticError, EPSGError,
              ETMError, FrechetError, GARSError, GeohashError, GeoidError, HausdorffError,
              HeightError, LazyImportError, LCCError, MGRSError, OSGRError, PGMError,
              PointsError, RefFrameError, SciPyError, SciPyWarning, TRFError, UnitError,
              UPSError, UTMError, UTMUPSError, VectorError, VincentyError,
              WebMercatorError, WGRSError):
        t.testError(E)
    t.testErrors(errors._InvalidError, errors._IsnotError)
    t.testKwds(errors._xkwds_get)
    t.testKwds(errors._xkwds_pop)
    t.results()
    t.exit()
