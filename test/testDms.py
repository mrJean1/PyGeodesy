
# -*- coding: utf-8 -*-

# Test degrees, minutes, seconds functions.

__all__ = ('Tests',)
__version__ = '20.05.12'

from base import TestsBase

from pygeodesy import F_D,   F_DM,   F_DMS,   F_DEG,   F_MIN,   F_SEC,   F_RAD, \
                      F_D_,  F_DM_,  F_DMS_,  F_DEG_,  F_MIN_,  F_SEC_,  F_RAD_, \
                      F_D__, F_DM__, F_DMS__, F_DEG__, F_MIN__, F_SEC__, F_RAD__, \
                      compassPoint, degDMS, fstr, parseDDDMMSS, parseDMS, \
                      ParseError, parse3llh, RangeError, rangerrors, toDMS


class Tests(TestsBase):

    def testDms(self):  # MCCABE 14
        for t in ('0.0°',
                  '0°',
              '''000°00'00"''',
              '''000°00'00.0"''',
              '''000° 00'00"''',
              '''000°00 ' 00.0"''',
              '''000° 00' 00.0'''):
            self.test("parseDMS('%s')" % (t,), parseDMS(t), '0.0')
        t = '000°-00′-00.0"'
        self.test("parseDMS('%s')" % (t,), parseDMS(t, sep='-'), '0.0')

        _X = ''
        for t, X, x in ((1, '1.0', _X),
                        (12, '12.0', _X),
                        (123, '123.0', _X),
                        (1234, '12.567', '1234.0'),
                        (12345, '123.75', '12345.0'),
                        (123456, '12.582', '123456.0'),
                        (1234567, '123.769', '1234567.0'),
                        (12345678, '1234.955', '12345678.0'),
                        (.1, 0.1, _X),
                        (1.2, 1.2, _X),
                        (12.3, 12.3, _X),
                        (123.4, 123.4, _X),
                        (1234.5, 12.575, 1234.5),
                        (12345.6, 123.76, 12345.6),
                        (123456.7, 12.582, 123456.7),
                        ('1N', '1.0', _X),
                        ('12S', '-12.0', _X),
                        ('012.3W', '-12.3', _X),
                        ('123E', '123.0', _X),
                        ('1234N', '12.567', '1234.0'),
                        ('12345E', '123.75', '12345.0'),
                        ('1234.5S', '-12.575', '-1234.5'),
                        ('12345.6E', '123.76', '12345.6'),
                        ('123456.7S', '-12.582', '-123456.7'),
                        ('1234567.8W', '-123.769', '-1234567.8'),
                        ('12345678E', '12345678.0', _X)):
            self.test('parseDDDMMSS(%r)' % (t,), fstr(parseDDDMMSS(t), prec=3), X)
            self.test('parseDMS(%r)' % (t,), fstr(parseDMS(t), prec=3), x or X)

        t = 345.0
        self.test('parseDDDMMSS(%r, NS)' % (t,), fstr(parseDDDMMSS(t, suffix='NS'), prec=3), '3.75')
        self.test('parseDDDMMSS(%r, EW)' % (t,), fstr(parseDDDMMSS(t, suffix='EW'), prec=3), '345.0')
        t = 5430.0
        self.test('parseDDDMMSS(%r, NS)' % (t,), fstr(parseDDDMMSS(t, suffix='NS'), prec=3), '54.5')
        self.test('parseDDDMMSS(%r, EW)' % (t,), fstr(parseDDDMMSS(t, suffix='EW'), prec=3), '54.5')

        for t, x in (('12E', '12.0'),
                     ('012.3S', '-12.3'),
                     ('123N', '123.0'),
                     ('1234E', '1234.0'),
                     ('12345N', '12345.0'),
                     ('1234.5W', '-1234.5'),
                     ('123456E', '123456.0'),
                     ('1234567S', '-1234567.0')):
            try:
                self.test('parseDDDMMSS(%r)' % (t,), parseDDDMMSS(t), ParseError.__name__)
            except ParseError as X:
                self.test('parseDDDMMSS(%r)' % (t,), repr(X), repr(X))
            self.test('parseDMS(%r)' % (t,), fstr(parseDMS(t), prec=5), x)

        r = rangerrors(True)
        try:
            self.test('parseDMS', parseDMS(181, clip=180), RangeError.__name__)
        except RangeError as x:
            self.test('parseDMS', str(x), str(x))
        rangerrors(False)
        try:
            self.test('parseDMS', parseDMS(-91, clip=90), RangeError.__name__)
        except RangeError as x:
            self.test('parseDMS', str(x), str(x))
        rangerrors(r)

        x = parse3llh('000° 00′ 05.31″W, 51° 28′ 40.12″ N')
        self.test('parse3llh', fstr(x, prec=6), '51.477811, -0.001475, 0.0')

        t = 'toDMS(%s)' % ('',)
        self.test(t, toDMS(45.99999, F_DM,  prec=1), '46°00.0′')  # not 45°60.0′
        self.test(t, toDMS(45.99999, F_DM,  prec=2), '46°00.0′')
        self.test(t, toDMS(45.9999,  F_DM,  prec=2), '45°59.99′')
        self.test(t, toDMS(45.99999, F_DM,  prec=3), '45°59.999′')
        self.test(t, toDMS(45.99999, F_DMS, prec=1), '46°00′00.0″')
        self.test(t, toDMS(45.99999, F_DMS, prec=2), '45°59′59.96″')
        self.test(t, toDMS(45.99999, F_DMS, prec=3), '45°59′59.964″')

        self.test(t, toDMS(45.76260),   '45°45′45.36″')

        for F, p, x in ((F_D,   None,   '45.7626°'),
                        (F_DM,  None,   "45°45.756'"),
                        (F_DMS, None, '''45°45'45.36"'''),
                        (F_DEG, None,   '45.7626'),
                        (F_MIN, None,   '4545.756'),
                        (F_SEC, None,   '454545.36'),
                        (F_RAD, None,   '0.79871'),

                        (F_D,    6,   '45.7626°'),
                        (F_DM,  -4,   "45°45.7560'"),
                        (F_DMS,  2, '''45°45'45.36"'''),
                        (F_DEG, -6,   '45.762600'),
                        (F_MIN, -5,   '4545.75600'),
                        (F_SEC, -3,   '454545.360'),
                        (F_RAD, -6,   '0.798708'),

                        (F_D_,    6,   '45.7626°'),
                        (F_DM_,  -4,   "45°45.7560'"),
                        (F_DMS_,  2, '''45°45'45.36"'''),
                        (F_DEG_, -6,   '45.762600'),
                        (F_MIN_, -5,   '4545.75600'),
                        (F_SEC_, -3,   '454545.360'),
                        (F_RAD_, -6,   '0.798708')):
            t = 'toDMS(%s)' % (F,)
            self.test(t, toDMS( 45.76260, F, prec=p),       x)
            self.test(t, toDMS(-45.76260, F, prec=p), '-' + x)

        for F, p, x in ((F_D__,    6,   '45.7626°'),
                        (F_DM__,  -4,   "45°45.7560'"),
                        (F_DMS__,  2, '''45°45'45.36"'''),
                        (F_DEG__, -6,   '45.762600'),
                        (F_MIN__, -5,   '4545.75600'),
                        (F_SEC__, -3,   '454545.360'),
                        (F_RAD__, -6,   '0.798708')):
            t = 'toDMS(%s)' % (F,)
            self.test(t, toDMS( 45.76260, F, prec=p), '+' + x)
            self.test(t, toDMS(-45.76260, F, prec=p), '-' + x)

        # <https://GitHub.com/chrisveness/geodesy/blob/master/test/dms-tests.js>
        for a, x in (((1,),   'N'),
                     ((0,),   'N'),
                     ((-1,),  'N'),
                     ((359,), 'N'),
                     ((24,),   'NNE'),
                     ((24, 1), 'N'),
                     ((24, 2), 'NE'),
                     ((24, 3), 'NNE'),
                     ((226,),   'SW'),
                     ((226, 1), 'W'),
                     ((226, 2), 'SW'),
                     ((226, 3), 'SW'),
                     ((237,),   'WSW'),
                     ((237, 1), 'W'),
                     ((237, 2), 'SW'),
                     ((237, 3), 'WSW'),
                     # Martin Schultz
                     ((11.25,), 'NNE'),
                     ((11.249,), 'N'),
                     ((-11.25,), 'N'),
                     ((348.749,), 'NNW'),
                     ((45, 1), 'E'),
                     ((44.99, 1), 'N'),
                     ((45, 2), 'NE'),
                     ((44.99, 2), 'NE'),
                     ((45, 3), 'NE'),
                     ((44.99, 3), 'NE'),
                     ((45, 4), 'NE'),
                     ((44.99, 4), 'NE'),  # XXX
                     ((22.5, 1), 'N'),
                     ((22.49, 1), 'N'),
                     ((22.5, 2), 'NE'),
                     ((22.49, 2), 'N'),
                     ((22.5, 3), 'NNE'),
                     ((22.49, 3), 'NNE'),
                     ((22.5, 4), 'NNE'),
                     ((22.49, 4), 'NNE'),  # XXX
                     ((11.25, 1), 'N'),
                     ((11.249, 1), 'N'),
                     ((11.25, 2), 'N'),
                     ((11.249, 2), 'N'),
                     ((11.25, 3), 'NNE'),
                     ((11.249, 3), 'N'),
                     ((11.25, 4), 'NbE'),
                     ((11.249, 4), 'NbE'),  # XXX
                     # examples
                     ((24, 1), 'N'),
                     ((24, 2), 'NE'),
                     ((24, 3), 'NNE'),
                     ((24,),   'NNE'),
                     ((18, 3), 'NNE'),
                     ((11, 4), 'NbE'),
                     ((30, 4), 'NEbN')):
            self.test('compassPoint%s' % (a,), compassPoint(*a), x)

        for a, x in enumerate(('NbE', 'NEbN', 'NEbE', 'EbN',
                               'EbS', 'SEbE', 'SEbS', 'SbE',
                               'SbW', 'SWbS', 'SWbW', 'WbS',
                               'WbN', 'NWbW', 'NWbN', 'NbW')):
            a = 11.25 + a * 22.5
            self.test('compassPoint(%s)' % (a,), compassPoint(a, 4), x)

        t = degDMS(1.0101, prec=5, s_D='', pos='+')  # coverage
        self.test('_DEG', t, '+1.0101')
        t = degDMS(0.0101, prec=5, s_S='', pos='+')
        self.test('_MIN', t, "+0.606'")
        t = degDMS(0.0101, prec=5, s_M='', pos='+')
        self.test('_SEC', t, '+36.36"')


if __name__ == '__main__':

    from pygeodesy import dms  # private

    t = Tests(__file__, __version__, dms)
    t.testDms()
    t.results()
    t.exit()
