
# -*- coding: utf-8 -*-

# Test degrees, minutes, seconds functions.

__all__ = ('Tests',)
__version__ = '18.09.28'

from base import TestsBase

from pygeodesy import F_D, F_DM, F_DMS, F_DEG, F_MIN, F_SEC, F_RAD, \
                      compassPoint, parse3llh, parseDMS, rangerrors, \
                      toDMS


class Tests(TestsBase):

    def testDms(self):
        # dms module tests
        self.test('parseDMS', parseDMS(  '0.0°'), '0.0')
        self.test('parseDMS', parseDMS(    '0°'), '0.0')
        self.test('parseDMS', parseDMS('''000°00'00"'''),   '0.0')
        self.test('parseDMS', parseDMS('''000°00'00.0"'''), '0.0')
        self.test('parseDMS', parseDMS('''000° 00'00"'''),    '0.0')
        self.test('parseDMS', parseDMS('''000°00 ' 00.0"'''), '0.0')

        r = rangerrors(True)
        try:
            self.test('parseDMS', parseDMS(181, clip=180), 'ValueError')
        except ValueError as x:
            self.test('parseDMS', str(x), '181.0 beyond 180 degrees')
        rangerrors(False)
        try:
            self.test('parseDMS', parseDMS(-91, clip=90), -90)
        except ValueError as x:
            self.test('parseDMS', str(x), '-90')
        rangerrors(r)

        x = parse3llh('000° 00′ 05.31″W, 51° 28′ 40.12″ N')
        x = ', '.join('%.6f' % a for a in x)  # XXX fStr
        self.test('parse3llh', x, '51.477811, -0.001475, 0.000000')

        for a, x in (((),            '''45°45'45.36"'''),
                     ((F_D,   None),   '45.7626°'),
                     ((F_DM,  None),   "45°45.756'"),
                     ((F_DMS, None), '''45°45'45.36"'''),
                     ((F_DEG, None),   '45.7626'),
                     ((F_MIN, None),   '4545.756'),
                     ((F_SEC, None),   '454545.36'),
                     ((F_RAD, None),   '0.79871'),
                     ((F_D,    6),   '45.7626°'),
                     ((F_DM,  -4),   "45°45.7560'"),
                     ((F_DMS,  2), '''45°45'45.36"'''),
                     ((F_DEG, -6),   '45.762600'),
                     ((F_MIN, -5),   '4545.75600'),
                     ((F_SEC, -3),   '454545.360'),
                     ((F_RAD, -6),   '0.798708')):
            t = 'toDMS(%s)' % (a[:1] or '')
            self.test(t, toDMS(45.76260, *a), x)

        # <http://GitHub.com/chrisveness/geodesy/blob/master/test/dms-tests.js>
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


if __name__ == '__main__':

    from pygeodesy import dms  # private

    t = Tests(__file__, __version__, dms)
    t.testDms()
    t.results()
    t.exit()
