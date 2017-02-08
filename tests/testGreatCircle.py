
# -*- coding: utf-8 -*-

# Transcribed from Objective-C GreatCircleTests by Brian Lambert (C)
# 2016 Softwarenerd at <https://github.com/softwarenerd/GreatCircle>

# The MIT License (MIT)
#
# Copyright (c) 2016 Softwarenerd.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# GreatCircleTests.m
# GreatCircleTests
#
# Created by Brian Lambert on 6/5/16.
# Copyright © 2016 Softwarenerd.

__version__ = '16.12.07'

if __name__ == '__main__':

    from tests import Tests as _Tests
    from geodesy import F_D, F_DMS, bearingDMS, \
                        sphericalNvector, sphericalTrigonometry

    class Tests(_Tests):
        # overload test() method
        def test(self, name, value, expect, fmt='%s', known=False):
            if isinstance(expect, float):
                expect = fmt % (expect,)
            _Tests.test(self, name, value, expect, fmt=fmt, known=known)

    def testGreatCircle(S):
        # run tests for spherical module S
        t = Tests(__file__ , __version__, S)

        # Indian Pond, in Piermond, NH.  My old Boy Scout Camp
        IndianPond = S.LatLon(43.930912, -72.053811)
        Eiffel     = S.LatLon(48.858158, 2.294825)
        Versailles = S.LatLon(48.804766, 2.120339)
        StGermain  = S.LatLon(48.897728, 2.094977)
        Orly       = S.LatLon(48.747114, 2.400526)

        # distance between the Eiffel Tower and Versailles in meter
        dEiffelToVersailles = 14084.280704919687
        mEiffelToVersailles = '%.4f'  # '%.8f'

        # initial and final bearings between Eiffel tower and Versailles
        ibEiffelToVersailles = 245.13460296861962
        fbEiffelToVersailles = 245.00325395138532
        fmtEiffelToVersailles = '%.8f'  # '%.14f'

        # initial and final bearing between Versailles and Eiffel tower
        ibVersaillesToEiffel = 65.003253951385318
        fbVersaillesToEiffel = 65.134602968619618
        fmtVersaillesToEiffel = '%.9f'  # '%.15f'

        # initial bearing for two locations that are the same
        b = IndianPond.bearingTo(IndianPond)
        t.test('InitialBearingSameLocations', b, 0.0, '%.1f')

        # initial bearing for two locations that are the equal
        b = IndianPond.bearingTo(IndianPond.copy())
        t.test('InitialBearingEqualLocations', b, 0.0, '%.1f')

        # final bearing for two locations that are the same
        b = IndianPond.finalBearingTo(IndianPond)
        t.test('FinalBearingSameLocations', b, 180.0, '%.1f')  # 0.0

        # final bearing for two locations that are the equal
        b = IndianPond.finalBearingTo(IndianPond.copy())
        t.test('FinalBearingEqualLocations', b, 180.0, '%.1f')  # 0.0

        # distance for two locations that are the same
        d = IndianPond.distanceTo(IndianPond)
        t.test('DistanceSameLocations', d, 0.0, '%.1f')

        # distance for two locations that are equal
        d = IndianPond.distanceTo(IndianPond.copy())
        t.test('DistanceEqualLocations', d, 0.0, '%.1f')

        # distance between Eiffel Tower and Versailles
        d = Eiffel.distanceTo(Versailles)
        t.test('DistanceEiffelToVersailles', d, dEiffelToVersailles, mEiffelToVersailles, known=True)

        # distance between Versailles and Eiffel Tower
        d = Versailles.distanceTo(Eiffel)
        t.test('DistanceVersaillesToEiffel', d, dEiffelToVersailles, mEiffelToVersailles, known=True)

        # initial bearing between Eiffel Tower and Versailles
        b = Eiffel.bearingTo(Versailles)
        t.test('InitialBearingEiffelToVersailles', b, ibEiffelToVersailles, fmtEiffelToVersailles)
        t.test('InitialBearingEiffelToVersailles(DMS)', bearingDMS(b, F_DMS, prec=4), '245°08′04.5707″')

        # initial bearing between Versailles and Eiffel Tower
        b = Versailles.bearingTo(Eiffel)
        t.test('InitialBearingVersaillesToEiffel', b, ibVersaillesToEiffel, fmtVersaillesToEiffel)
        t.test('InitialBearingVersaillesToEiffel(DMS)', bearingDMS(b, F_DMS, prec=4), '065°00′11.7142″')

        # final bearing between Eiffel Tower and Versailles
        b = Eiffel.finalBearingTo(Versailles)
        t.test('FinalBearingEiffelToVersailles', b, fbEiffelToVersailles, fmtEiffelToVersailles)
        t.test('FinalBearingEiffelToVersailles(DMS)', bearingDMS(b, F_DMS, prec=4), '245°00′11.7142″')

        # final bearing between Versailles and Eiffel Tower
        b = Versailles.finalBearingTo(Eiffel)
        t.test('FinalBearingVersaillesToEiffel', b, fbVersaillesToEiffel, fmtVersaillesToEiffel)
        t.test('FinalBearingVersaillesToEiffel(DMS)', bearingDMS(b, F_DMS, prec=4), '065°08′04.5707″')

        # generating a location for Versailles based on bearing and distance
        v = Eiffel.destination(dEiffelToVersailles, ibEiffelToVersailles)
        t.test('GenerateLocationVersailles', v, str(Versailles))

        # generating a location for Eiffel based on bearing and distance.
        e = Versailles.destination(dEiffelToVersailles, ibVersaillesToEiffel)
        t.test('GenerateLocationEiffel', e, str(Eiffel))

        # midpoint between the Eiffel and Versailles
        a = Eiffel.midpointTo(Versailles)
        b = Eiffel.destination(dEiffelToVersailles / 2.0, ibEiffelToVersailles)
        t.test('MidpointEiffelToVersailles', a, str(b))
        t.test('MidpointEiffelToVersailles(DMS)', a.toStr(F_DMS, prec=4), '48°49′53.3817″N, 002°12′27.1279″E')
        a = Eiffel.distanceTo(a)
        b = Versailles.distanceTo(b)
        t.test('MidpointEiffelToVersailles(m)', a, str(b), known=True)

        # midpoint between Versailles and the Eiffel Tower
        a = Versailles.midpointTo(Eiffel)
        b = Versailles.destination(dEiffelToVersailles / 2.0,
                                   ibVersaillesToEiffel)
        t.test('MidpointVersaillesToEiffel', a, str(b), known=True)
        t.test('MidpointVersaillesToEiffel(DMS)', a.toStr(F_DMS, prec=4), '48°49′53.3817″N, 002°12′27.1279″E')
        a = Versailles.distanceTo(a)
        b = Eiffel.distanceTo(b)
        t.test('MidpointVersaillesToEiffel(m)', a, str(b), known=True)

        # intersection.
        b = StGermain.bearingTo(Orly)
        i = StGermain.intersection(b, Eiffel, ibEiffelToVersailles)
        t.test('Intersection', i.toStr(F_D, prec=9), '48.83569095°N, 002.221252031°E')  # '48.83569094988361°N, ...
        t.test('Intersection', i.toStr(F_D, prec=13), '48.8356909498836°N, 002.2212520313074°E')  # 002.2212520313073583°E'

        # cross-track distance test of a point 90° and 200 meters away
        m = Eiffel.midpointTo(Versailles)
        b = Eiffel.bearingTo(Versailles)
        p = m.destination(200.0, (b + 90) % 360.0)
        d = p.crossTrackDistanceTo(Eiffel, Versailles)
        t.test('CrossTrackDistance90Degrees200Meters', d, 200.0, '%0.1f')

        # cross-track distance test of a point 270° and 200 meters away
        m = Eiffel.midpointTo(Versailles)
        b = Eiffel.bearingTo(Versailles)
        p = m.destination(200.0, (b + 270) % 360.0)
        d = p.crossTrackDistanceTo(Eiffel, Versailles)
        t.test('CrossTrackDistance270Degrees200Meters', d, -200.0, '%0.1f')

        # cross-track distance that should be very close to 0
        m = Eiffel.midpointTo(Versailles)
        d = abs(m.crossTrackDistanceTo(Eiffel, Versailles))
        t.test('CrossTrackDistanceThatShouldBeVeryCloseToZero', d, '0.00000000', '%.8f')

        t.results()
        return t

    e = testGreatCircle(sphericalNvector).errors()
    t = testGreatCircle(sphericalTrigonometry)
    t.exit(e)

    # Typical test results (on MacOS 10.12.2):

    # testing geodesy.sphericalNvector version 17.02.07
    # test 1 InitialBearingSameLocations: 0.0
    # test 2 InitialBearingEqualLocations: 0.0
    # test 3 FinalBearingSameLocations: 180.0
    # test 4 FinalBearingEqualLocations: 180.0
    # test 5 DistanceSameLocations: 0.0
    # test 6 DistanceEqualLocations: 0.0
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 9 InitialBearingEiffelToVersailles: 245.13460297
    # test 10 InitialBearingEiffelToVersailles(DMS): 245°08′04.5707″
    # test 11 InitialBearingVersaillesToEiffel: 65.003253951
    # test 12 InitialBearingVersaillesToEiffel(DMS): 065°00′11.7142″
    # test 13 FinalBearingEiffelToVersailles: 245.00325395
    # test 14 FinalBearingEiffelToVersailles(DMS): 245°00′11.7142″
    # test 15 FinalBearingVersaillesToEiffel: 65.134602969
    # test 16 FinalBearingVersaillesToEiffel(DMS): 065°08′04.5707″
    # test 17 GenerateLocationVersailles: 48.804766°N, 002.120339°E
    # test 18 GenerateLocationEiffel: 48.858158°N, 002.294825°E
    # test 19 MidpointEiffelToVersailles: 48.831495°N, 002.207536°E
    # test 20 MidpointEiffelToVersailles(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 21 MidpointEiffelToVersailles(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 23 MidpointVersaillesToEiffel(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 24 MidpointVersaillesToEiffel(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 25 Intersection: 48.83569095°N, 002.221252031°E
    # test 26 Intersection: 48.8356909498836°N, 002.2212520313074°E
    # test 27 CrossTrackDistance90Degrees200Meters: 200.0
    # test 28 CrossTrackDistance270Degrees200Meters: -200.0
    # test 29 CrossTrackDistanceThatShouldBeVeryCloseToZero: 0.00000000
    # 5 geodesy.sphericalNvector tests (17.2%) FAILED, incl. 5 KNOWN (Python 2.7.13 64bit)

    # testing geodesy.sphericalTrigonometry version 17.02.07
    # test 1 InitialBearingSameLocations: 0.0
    # test 2 InitialBearingEqualLocations: 0.0
    # test 3 FinalBearingSameLocations: 180.0
    # test 4 FinalBearingEqualLocations: 180.0
    # test 5 DistanceSameLocations: 0.0
    # test 6 DistanceEqualLocations: 0.0
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 9 InitialBearingEiffelToVersailles: 245.13460297
    # test 10 InitialBearingEiffelToVersailles(DMS): 245°08′04.5707″
    # test 11 InitialBearingVersaillesToEiffel: 65.003253951
    # test 12 InitialBearingVersaillesToEiffel(DMS): 065°00′11.7142″
    # test 13 FinalBearingEiffelToVersailles: 245.00325395
    # test 14 FinalBearingEiffelToVersailles(DMS): 245°00′11.7142″
    # test 15 FinalBearingVersaillesToEiffel: 65.134602969
    # test 16 FinalBearingVersaillesToEiffel(DMS): 065°08′04.5707″
    # test 17 GenerateLocationVersailles: 48.804766°N, 002.120339°E
    # test 18 GenerateLocationEiffel: 48.858158°N, 002.294825°E
    # test 19 MidpointEiffelToVersailles: 48.831495°N, 002.207536°E
    # test 20 MidpointEiffelToVersailles(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 21 MidpointEiffelToVersailles(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 23 MidpointVersaillesToEiffel(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 24 MidpointVersaillesToEiffel(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 25 Intersection: 48.83569095°N, 002.221252031°E
    # test 26 Intersection: 48.8356909498836°N, 002.2212520313074°E
    # test 27 CrossTrackDistance90Degrees200Meters: 200.0
    # test 28 CrossTrackDistance270Degrees200Meters: -200.0
    # test 29 CrossTrackDistanceThatShouldBeVeryCloseToZero: 0.00000000
    # 5 geodesy.sphericalTrigonometry tests (17.2%) FAILED, incl. 5 KNOWN (Python 2.7.13 64bit)

    # testing sphericalNvector version 17.02.07
    # test 1 InitialBearingSameLocations: 0.0
    # test 2 InitialBearingEqualLocations: 0.0
    # test 3 FinalBearingSameLocations: 180.0
    # test 4 FinalBearingEqualLocations: 180.0
    # test 5 DistanceSameLocations: 0.0
    # test 6 DistanceEqualLocations: 0.0
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 9 InitialBearingEiffelToVersailles: 245.13460297
    # test 10 InitialBearingEiffelToVersailles(DMS): 245°08′04.5707″
    # test 11 InitialBearingVersaillesToEiffel: 65.003253951
    # test 12 InitialBearingVersaillesToEiffel(DMS): 065°00′11.7142″
    # test 13 FinalBearingEiffelToVersailles: 245.00325395
    # test 14 FinalBearingEiffelToVersailles(DMS): 245°00′11.7142″
    # test 15 FinalBearingVersaillesToEiffel: 65.134602969
    # test 16 FinalBearingVersaillesToEiffel(DMS): 065°08′04.5707″
    # test 17 GenerateLocationVersailles: 48.804766°N, 002.120339°E
    # test 18 GenerateLocationEiffel: 48.858158°N, 002.294825°E
    # test 19 MidpointEiffelToVersailles: 48.831495°N, 002.207536°E
    # test 20 MidpointEiffelToVersailles(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 21 MidpointEiffelToVersailles(m): 7042.150047881914  FAILED, KNOWN, expected 7042.159743304898
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 23 MidpointVersaillesToEiffel(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 24 MidpointVersaillesToEiffel(m): 7042.150047882616  FAILED, KNOWN, expected 7042.159743304495
    # test 25 Intersection: 48.83569095°N, 002.221252031°E
    # test 26 Intersection: 48.8356909498836°N, 002.2212520313074°E
    # test 27 CrossTrackDistance90Degrees200Meters: 200.0
    # test 28 CrossTrackDistance270Degrees200Meters: -200.0
    # test 29 CrossTrackDistanceThatShouldBeVeryCloseToZero: 0.00000000
    # 5 sphericalNvector tests (17.2%) FAILED, incl. 5 KNOWN (Python 3.6.0 64bit)

    # testing sphericalTrigonometry version 17.02.07
    # test 1 InitialBearingSameLocations: 0.0
    # test 2 InitialBearingEqualLocations: 0.0
    # test 3 FinalBearingSameLocations: 180.0
    # test 4 FinalBearingEqualLocations: 180.0
    # test 5 DistanceSameLocations: 0.0
    # test 6 DistanceEqualLocations: 0.0
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 9 InitialBearingEiffelToVersailles: 245.13460297
    # test 10 InitialBearingEiffelToVersailles(DMS): 245°08′04.5707″
    # test 11 InitialBearingVersaillesToEiffel: 65.003253951
    # test 12 InitialBearingVersaillesToEiffel(DMS): 065°00′11.7142″
    # test 13 FinalBearingEiffelToVersailles: 245.00325395
    # test 14 FinalBearingEiffelToVersailles(DMS): 245°00′11.7142″
    # test 15 FinalBearingVersaillesToEiffel: 65.134602969
    # test 16 FinalBearingVersaillesToEiffel(DMS): 065°08′04.5707″
    # test 17 GenerateLocationVersailles: 48.804766°N, 002.120339°E
    # test 18 GenerateLocationEiffel: 48.858158°N, 002.294825°E
    # test 19 MidpointEiffelToVersailles: 48.831495°N, 002.207536°E
    # test 20 MidpointEiffelToVersailles(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 21 MidpointEiffelToVersailles(m): 7042.150047881641  FAILED, KNOWN, expected 7042.159743304608
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 23 MidpointVersaillesToEiffel(DMS): 48°49′53.3817″N, 002°12′27.1279″E
    # test 24 MidpointVersaillesToEiffel(m): 7042.150047882363  FAILED, KNOWN, expected 7042.159743304893
    # test 25 Intersection: 48.83569095°N, 002.221252031°E
    # test 26 Intersection: 48.8356909498836°N, 002.2212520313074°E
    # test 27 CrossTrackDistance90Degrees200Meters: 200.0
    # test 28 CrossTrackDistance270Degrees200Meters: -200.0
    # test 29 CrossTrackDistanceThatShouldBeVeryCloseToZero: 0.00000000
    # 5 sphericalTrigonometry tests (17.2%) FAILED, incl. 5 KNOWN (Python 3.6.0 64bit)
