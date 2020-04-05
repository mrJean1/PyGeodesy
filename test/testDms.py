
# -*- coding: utf-8 -*-

# Test degrees, minutes, seconds functions.

__all__ = ('Tests',)
__version__ = '20.04.04'

from base import TestsBase

from pygeodesy import F_D,   F_DM,   F_DMS,   F_DEG,   F_MIN,   F_SEC,   F_RAD, \
                      F_D_,  F_DM_,  F_DMS_,  F_DEG_,  F_MIN_,  F_SEC_,  F_RAD_, \
                      F_D__, F_DM__, F_DMS__, F_DEG__, F_MIN__, F_SEC__, F_RAD__, \
                      compassPoint, degDMS, fStr, parse3llh, parseDMS, rangerrors, \
                      toDMS


class Tests(TestsBase):

    def testDms(self):
        # dms module tests
        for i, t in enumerate(('0.0°',
                               '0°',
                           '''000°00'00"''',
                           '''000°00'00.0"''',
                           '''000° 00'00"''',
                           '''000°00 ' 00.0"''',
                           '''000° 00' 00.0''')):
            self.test('parseDMS' + str(i + 1), parseDMS(t), '0.0')
        self.test('parseDMS' + str(i + 2), parseDMS('000°-00′-00.0"', sep='-'), '0.0')

        r = rangerrors(True)
        try:
            self.test('parseDMS', parseDMS(181, clip=180), 'ValueError')
        except ValueError as x:
            self.test('parseDMS', str(x), '181.0 beyond 180 degrees')
        rangerrors(False)
        try:
            self.test('parseDMS', parseDMS(-91, clip=90), '-90.0')
        except ValueError as x:
            self.test('parseDMS', str(x), '-90.0')
        rangerrors(r)

        x = parse3llh('000° 00′ 05.31″W, 51° 28′ 40.12″ N')
        self.test('parse3llh', fStr(x, prec=6), '51.477811, -0.001475, 0.0')

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
