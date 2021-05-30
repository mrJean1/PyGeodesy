
# -*- coding: utf-8 -*-

# Transcoded from Objective-C GreatCircleTests by Brian Lambert (C)
# 2016 Softwarenerd at <https://GitHub.com/softwarenerd/GreatCircle>

# The MIT License (MIT)
#
# Copyright (C) 2016 Softwarenerd.
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
# Copyright (C) 2016 Softwarenerd.

__all__ = ()
__version__ = '20.09.09'  # '18.10.02'

from base import TestsBase

from pygeodesy import F_D, F_DMS, classname, bearingDMS, CrossError, crosserrors


class Tests(TestsBase):

    def testGreatCircle(self, module):  # spherical only

        self.subtitle(module, 'GreatCircle')

        LatLon = module.LatLon

        # Indian Pond, in Piermond, NH.  My old Boy Scout Camp
        IndianPond = LatLon(43.930912, -72.053811)
        Eiffel     = LatLon(48.858158, 2.294825)
        Versailles = LatLon(48.804766, 2.120339)
        StGermain  = LatLon(48.897728, 2.094977)
        Orly       = LatLon(48.747114, 2.400526)

        # distance between the Eiffel Tower and Versailles in meter
        dEiffelToVersailles = 14084.280704919687
        xEiffelToVersailles = '%.6f'  # '%.8f'

        # initial and final bearings between Eiffel tower and Versailles
        ibEiffelToVersailles = 245.13460296861962
        fbEiffelToVersailles = 245.00325395138532
        xbEiffelToVersailles = '%.8f'  # '%.14f'

        # initial and final bearing between Versailles and Eiffel tower
        ibVersaillesToEiffel = 65.003253951385318
        fbVersaillesToEiffel = 65.134602968619618
        xbVersaillesToEiffel = '%.9f'  # '%.15f'

        xMidpoint = '%.8f'

        c = crosserrors(False)  # no CrossErrors!

        # initial bearing for two locations that are the same
        b = IndianPond.initialBearingTo(IndianPond)
        self.test('InitialBearingSameLocations', b, 0.0, prec=1)

        # initial bearing for two locations that are the equal
        b = IndianPond.initialBearingTo(IndianPond.copy())
        self.test('InitialBearingEqualLocations', b, 0.0, prec=1)

        # final bearing for two locations that are the same
        b = IndianPond.finalBearingTo(IndianPond)
        self.test('FinalBearingSameLocations', b, 180.0, prec=1)  # 0.0

        # final bearing for two locations that are the equal
        b = IndianPond.finalBearingTo(IndianPond.copy())
        self.test('FinalBearingEqualLocations', b, 180.0, prec=1)  # 0.0

        c = crosserrors(c)  # with CrossErrors!

        try:  # should raise CrossError
            b = str(IndianPond.initialBearingTo(IndianPond, raiser=True))
        except CrossError as x:
            b = str(x)
        self.test('FinalBearingCrossError', b, 'points (%s(43°55′51.28″N, 072°03′13.72″W)): coincident' % (classname(IndianPond),))

        # distance for two locations that are the same
        d = IndianPond.distanceTo(IndianPond)
        self.test('DistanceSameLocations', d, 0.0, prec=1)

        # distance for two locations that are equal
        d = IndianPond.distanceTo(IndianPond.copy())
        self.test('DistanceEqualLocations', d, 0.0, prec=1)

        # distance between Eiffel Tower and Versailles
        d = Eiffel.distanceTo(Versailles)
        self.test('DistanceEiffelToVersailles', d, dEiffelToVersailles, fmt=xEiffelToVersailles, known=True)

        # distance between Versailles and Eiffel Tower
        d = Versailles.distanceTo(Eiffel)
        self.test('DistanceVersaillesToEiffel', d, dEiffelToVersailles, fmt=xEiffelToVersailles, known=True)

        # initial bearing between Eiffel Tower and Versailles
        b = Eiffel.initialBearingTo(Versailles)
        self.test('InitialBearingEiffelToVersailles', b, ibEiffelToVersailles, fmt=xbEiffelToVersailles)
        self.test('InitialBearingEiffelToVersailles(DMS)', bearingDMS(b, F_DMS, prec=4), '245°08′04.5707″')

        # initial bearing between Versailles and Eiffel Tower
        b = Versailles.initialBearingTo(Eiffel)
        self.test('InitialBearingVersaillesToEiffel', b, ibVersaillesToEiffel, fmt=xbVersaillesToEiffel)
        self.test('InitialBearingVersaillesToEiffel(DMS)', bearingDMS(b, F_DMS, prec=4), '65°00′11.7142″')

        # final bearing between Eiffel Tower and Versailles
        b = Eiffel.finalBearingTo(Versailles)
        self.test('FinalBearingEiffelToVersailles', b, fbEiffelToVersailles, fmt=xbEiffelToVersailles)
        self.test('FinalBearingEiffelToVersailles(DMS)', bearingDMS(b, F_DMS, prec=4), '245°00′11.7142″')

        # final bearing between Versailles and Eiffel Tower
        b = Versailles.finalBearingTo(Eiffel)
        self.test('FinalBearingVersaillesToEiffel', b, fbVersaillesToEiffel, fmt=xbVersaillesToEiffel)
        self.test('FinalBearingVersaillesToEiffel(DMS)', bearingDMS(b, F_DMS, prec=4), '65°08′04.5707″')

        # generating a location for Versailles based on bearing and distance
        v = Eiffel.destination(dEiffelToVersailles, ibEiffelToVersailles)
        self.test('GenerateLocationVersailles', v, str(Versailles))

        # generating a location for Eiffel based on bearing and distance.
        e = Versailles.destination(dEiffelToVersailles, ibVersaillesToEiffel)
        self.test('GenerateLocationEiffel', e, str(Eiffel))

        # midpoint between the Eiffel and Versailles
        a = Eiffel.midpointTo(Versailles)
        b = Eiffel.destination(dEiffelToVersailles / 2.0, ibEiffelToVersailles)
        self.test('MidpointEiffelToVersailles', a, str(b))
        self.test('MidpointEiffelToVersailles(DMS)', a.toStr(F_DMS, prec=4), '48°49′53.3817″N, 002°12′27.1279″E')
        a = Eiffel.distanceTo(a)
        m = Versailles.distanceTo(b)
        self.test('MidpointEiffelToVersailles(m)', a, m, fmt=xMidpoint, known=True)

        # midpoint between Versailles and the Eiffel Tower
        a = Versailles.midpointTo(Eiffel)
        b = Versailles.destination(dEiffelToVersailles / 2.0, ibVersaillesToEiffel)
        self.test('MidpointVersaillesToEiffel', a, str(b), known=True)
        self.test('MidpointVersaillesToEiffel(DMS)', a.toStr(F_DMS, prec=4), '48°49′53.3817″N, 002°12′27.1279″E')
        a = Versailles.distanceTo(a)
        m = Eiffel.distanceTo(b)
        self.test('MidpointVersaillesToEiffel(m)', a, m, fmt=xMidpoint, known=True)

        # intersection.
        b = StGermain.initialBearingTo(Orly)
        i = StGermain.intersection(b, Eiffel, ibEiffelToVersailles)
        self.test('Intersection', i.toStr(F_D, prec=9), '48.83569095°N, 002.221252031°E')  # '48.83569094988361°N, ...
        self.test('Intersection', i.toStr(F_D, prec=13), '48.8356909498836°N, 002.2212520313074°E')  # 002.2212520313073583°E'

        # cross-track distance test of a point 90° and 200 meters away
        m = Eiffel.midpointTo(Versailles)
        b = Eiffel.initialBearingTo(Versailles)
        p = m.destination(200.0, (b + 90) % 360.0)
        d = p.crossTrackDistanceTo(Eiffel, Versailles)
        self.test('CrossTrackDistance200m+90°', d, 200.0, prec=1)

        # cross-track distance test of a point 270° and 200 meters away
        m = Eiffel.midpointTo(Versailles)
        b = Eiffel.initialBearingTo(Versailles)
        p = m.destination(200.0, (b + 270) % 360.0)
        d = p.crossTrackDistanceTo(Eiffel, Versailles)
        self.test('CrossTrackDistance200m+270°', d, -200.0, prec=1)

        # cross-track distance that should be very close to 0
        m = Eiffel.midpointTo(Versailles)
        d = abs(m.crossTrackDistanceTo(Eiffel, Versailles))
        self.test('CrossTrackDistanceCloseToZero', d, '0.0000000', fmt='%.7f')


if __name__ == '__main__':

    from pygeodesy import sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)

    t.testGreatCircle(sphericalNvector)

    t.testGreatCircle(sphericalTrigonometry)

    t.results()
    t.exit()
