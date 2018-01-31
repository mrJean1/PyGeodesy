
# -*- coding: utf-8 -*-

# Test degrees, minutes, seconds functions.

__all__ = ('Tests',)
__version__ = '18.01.31'

from base import TestsBase

from pygeodesy import F_D, F_DM, F_DMS, F_DEG, F_MIN, F_SEC, F_RAD, \
                      compassAngle, compassPoint, equirectangular, \
                      m2km, parse3llh, parseDMS, rangerrors, toDMS

from pygeodesy.sphericalTrigonometry import LatLon
_LHR = LatLon(51.47,   0.4543)
_FRA = LatLon(50.0379, 8.5622)


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

        self.test('compassAngle0', compassAngle(0, 0,  0,  0),   0.0)
        self.test('compassAngle1', compassAngle(0, 0, 10, 10),  45.0)
        self.test('compassAngle2', compassAngle(0, 0,  0, 10),  90.0)
        self.test('compassAngle3', compassAngle(0, 0, -1,  0), 180.0)
        self.test('compassAngle4', compassAngle(0, 0, -1, -1), 225.0)
        self.test('compassAngle5', compassAngle(0, 0,  0, -1), 270.0)
        self.test('compassAngle6', compassAngle(0, 0,  1, -1), 315.0)
        self.test('compassAngle7', compassAngle(0, 0, 99, -1), 359.4, fmt='%.1f')

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
                     ((237, 3), 'WSW')):
            self.test('compassPoint', compassPoint(*a), x)

        for a, x in enumerate(('NbE', 'NEbN', 'NEbE', 'EbN',
                               'EbS', 'SEbE', 'SEbS', 'SbE',
                               'SbW', 'SWbS', 'SWbW', 'WbS',
                               'WbN', 'NWbW', 'NWbN', 'NbW')):
            a = 11.25 + a * 22.5
            self.test('compassPoint', compassPoint(a, 4), x)

        a = compassAngle(_LHR.lat, _LHR.lon, _FRA.lat, _FRA.lon)
        self.test('compassAngle', a, 100.016848, fmt='%.6f')
        b = _LHR.initialBearingTo(_FRA)
        self.test('initialBearingTo', b, 102.432182, fmt='%.6f')
        d = equirectangular(_LHR.lat, _LHR.lon, _FRA.lat, _FRA.lon)
        self.test('equirectangular', m2km(d), 592.185, fmt='%.3f')
        d = _LHR.distanceTo(_FRA)
        self.test('distanceTo', m2km(d), 591.831, fmt='%.3f')


if __name__ == '__main__':

    from pygeodesy import dms  # private

    t = Tests(__file__, __version__, dms)
    t.testDms()
    t.results(nl=0)
    t.exit()
