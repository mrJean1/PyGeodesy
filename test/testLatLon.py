
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '20.06.15'

from base import geographiclib, TestsBase

from pygeodesy import R_NM, F_D, F_DM, F_DMS, F_RAD, \
                      degrees, fStr, isclockwise, isconvex, \
                      isenclosedBy, ispolar, m2km, m2NM, \
                      VincentyError  # PYCHOK expected
from pygeodesy.named import Bounds2Tuple, \
                            LatLon2Tuple, LatLon3Tuple, \
                            PhiLam2Tuple, PhiLam3Tuple, \
                            Vector3Tuple, Vector4Tuple


class Tests(TestsBase):

    def testLatLon(self, module, Sph=False, Nv=True):  # MCCABE 45

        self.subtitle(module, 'LatLon')

        LatLon = module.LatLon

        # basic LatLon class tests
        p = LatLon(52.20472, 0.14056)
        self.test('isEllipsoidal', p.isEllipsoidal, not Sph)
        self.test('isSpherical', p.isSpherical, Sph)

        self.test('lat/lonDMS', p, '52.20472°N, 000.14056°E')  # 52.20472°N, 000.14056°E
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 3),  '''52°12.283'N, 000°08.434'E''')
        self.test('lat/lonDMS F_DM', p.toStr(F_DM, 4),  '''52°12.2832'N, 000°08.4336'E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 0), '''52°12'17"N, 000°08'26"E''')
        self.test('lat/lonDMS F_DMS', p.toStr(F_DMS, 1), '''52°12'17.0"N, 000°08'26.0"E''')
        self.test('lat/lonDMS F_RAD', p.toStr(F_RAD, 6), '0.911144N, 0.002453E')

        q = LatLon(*map(degrees, p.philam))
        self.test('isequalTo', q.isequalTo(p), True)
        q = LatLon(*map(degrees, p.to2ab()))
        self.test('isequalTo', q.isequalTo(p), True)

        self.test('latlon2', fStr(q.latlon2(5)), '52.20472, 0.14056')
        self.test('latlon2', fStr(q.latlon2(4)), '52.2047, 0.1406')
        self.test('latlon2', fStr(q.latlon2(3)), '52.205, 0.141')
        self.test('latlon2', fStr(q.latlon2(2)), '52.2, 0.14')
        self.test('latlon2', fStr(q.latlon2(1)), '52.2, 0.1')
        self.test('latlon2', fStr(q.latlon2()),  '52.0, 0.0')

        FRA = LatLon(50.0379, 8.5622)
        LHR = LatLon(51.47, 0.4543)
        # <https://www.EdWilliams.org/avform.htm#XTE>
        JFK = LatLon(degrees(0.709186), -degrees(1.287762))
        LAX = LatLon(33.+57./60, -(118.+24./60))
        Rav = m2NM(6366710)  # avg. earth radius in NM
        # <https://GeographicLib.SourceForge.io/html/python/examples.html>
        WNZ = LatLon(-41.32, 174.81, name='Wellington, NZ')
        SAL = LatLon(40.96, -5.50, name='Salamanca, ES')
        BJS = LatLon(40.1, 116.6)  # Beijing Airport
        SFO = LatLon(37.6, -122.4)  # San Francisco

        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        self.test('isequalTo', p.isequalTo(q), False)

        a = p.antipode()
        self.test('antipode1', a, '52.205°S, 179.881°W')
        self.test('antipode2', a.isantipodeTo(p), True)
        b = a.antipode()
        self.test('antipode3', b, '52.205°N, 000.119°E')
        self.test('antipode4', a.isantipodeTo(b), True)
        self.test('antipode5', b, p)

        if hasattr(LatLon, 'initialBearingTo'):
            b = p.initialBearingTo(q)
            self.test('initialBearingTo', b, 156.1666 if Sph else 156.1106, fmt='%.4f')  # 156.2
            b = p.finalBearingTo(q)
            self.test('finalBearingTo', b, 157.8904 if Sph else 157.8345, fmt='%.4f')
            b = LAX.initialBearingTo(JFK)
            self.test('initialBearingTo', b, 65.8921 if Sph else 65.9335, fmt='%.4f')  # PYCHOK test attr?
            b = LAX.finalBearingTo(JFK)
            self.test('finalBearingTo', b, 93.8581 if Sph else 93.9034, fmt='%.4f')  # PYCHOK test attr?

        if hasattr(LatLon, 'bearingTo2'):
            b = p.bearingTo2(q)
            self.test('bearingTo2', fStr(b, prec=4), '156.1666, 157.8904' if Sph else '156.1106, 157.8345')  # 156.2
            # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
            b = LatLon(85, -135), LatLon(85, -45), LatLon(85, 45), LatLon(85, 135), LatLon(85, -135)
            self.test('ispolar', ispolar(b), True)  # PYCHOK test attr?

        c = p.copy()
        self.test('copy', p.isequalTo(c), True)
        self.test('__eq__', p == c, True)
        self.test('__ne__', p != c, False)

        d = p.equirectangularTo(q)
        self.test('equirectangularTo', d, '404329.56', fmt='%.2f')

        if hasattr(LatLon, 'distanceTo'):
            d = p.distanceTo(q)
            self.test('distanceTo', d, '404279.720589' if Sph else ('404279.720589' if Nv else '404607.805988'), fmt='%.6f')  # 404300
            d = q.distanceTo(p)
            self.test('distanceTo', d, '404279.720589' if Sph else ('404279.720589' if Nv else '404607.805988'), fmt='%.6f')  # 404300
            d = LAX.distanceTo(JFK, radius=R_NM) if Sph else LAX.distanceTo(JFK)
            self.test('distanceTo', d, 2145 if Sph else (3972863 if Nv else 3981601), fmt='%.0f')  # PYCHOK test attr?
            if not Nv:  # <https://GeographicLib.SourceForge.io/html/python/examples.html>
                self.test('antipodal', WNZ.isantipodeTo(SAL, eps=0.1), False)
                try:
                    d = WNZ.distanceTo(SAL, wrap=False)
                    self.test('distanceTo dateline', d, 19119590.551 if Sph else 19959679.267, fmt='%.3f', known=True)  # PYCHOK test attr?
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo dateline', x, 'no convergence ...', known=True)
                try:
                    d = WNZ.distanceTo(SAL, wrap=True)
                    self.test('distanceTo unrolled', d, 19119590.551 if Sph else 19959679.267, fmt='%.3f', known=True)  # PYCHOK test attr?
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo unrolled', x, 'no convergence ...', known=True)
                self.test('antipodal', BJS.isantipodeTo(SFO, eps=0.1), False)
                d = BJS.distanceTo(SFO, wrap=False)
                self.test('distanceTo dateline', d, 9491735 if Sph else 9513998, fmt='%.0f')  # PYCHOK test attr?
                d = BJS.distanceTo(SFO, wrap=True)
                self.test('distanceTo unrolled', d, 9491735 if Sph else 9513998, fmt='%.0f')  # PYCHOK test attr?

            # <https://GitHub.com/chrisveness/geodesy/issues/64>
            d = LatLon(20, 0).distanceTo(LatLon(-2, 180))
            self.test('distanceTo', d, 18013602.92 if Sph or Nv else 18012714.66, fmt='%.2f')
            try:
                d = LatLon(0, 0).distanceTo(LatLon(0, 180))  # antipodal
                self.test('distanceTo', d, 20015114.35 if Sph else 20003931.46, fmt='%.2f', known=Nv)  # PYCHOK 0.0 for Nv
            except ValueError as x:
                x = str(x)
                self.test('distanceTo', x, 'ambiguous, antipodal ...', known=True)

        if hasattr(LatLon, 'distanceTo3'):  # and not Sph:
            for w in (False, True):
                try:
                    d = WNZ.distanceTo3(SAL, wrap=w)  # XXX expected values?
                    self.test('distanceTo3 dateline', fStr(d, prec=4), '19959679.2674, 161.0677, 18.8252', known=True)  # PYCHOK test attr?
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo3 dateline', x, 'no convergence ...', known=True)
                d = BJS.distanceTo3(SFO, wrap=w)
                self.test('distanceTo3 dateline', fStr(d, prec=4), '9513997.9901, 42.9164, 138.8903')  # PYCHOK test attr?

        if hasattr(LatLon, 'intermediateTo'):
            i = p.intermediateTo(q, 0.25)
            self.test('intermediateTo', i, '51.372084°N, 000.707337°E' if Sph
                                      else '51.372294°N, 000.707192°E')
            self.test('intermediateTo', isinstance(i, LatLon), True)

            if hasattr(p, 'distanceTo'):
                d = p.distanceTo(q)
                self.test('intermediateTo', d, '404279.721', fmt='%.3f')  # PYCHOK test attr?

            i = p.intermediateTo(q, 5)
            self.test('intermediateTo+5', i, '35.160975°N, 008.989542°E' if Sph
                                        else '35.560239°N, 008.833512°E')
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo+5', p.distanceTo(i) / d, '4.885' if Nv and not Sph else '5.000', fmt='%.3f')  # PYCHOK test attr?

            i = p.intermediateTo(q, -4)
            self.test('intermediateTo-4', i, '64.911647°N, 013.726301°W' if Sph
                                        else '64.570387°N, 013.156352°W')
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo-4', p.distanceTo(i) / d, '3.885' if Nv and not Sph else '4.000', fmt='%.3f')  # PYCHOK test attr?

            # courtesy of <https://GitHub.com/bakakaldsas>
            i = LatLon(52, 0, 100).intermediateTo(LatLon(48, 2, 200), 0.25)
            self.test('intermediateTo-h', i.height, '125.000', fmt='%.3f')  # PYCHOK test attr?

        if hasattr(LatLon, 'intermediateChordTo'):
            i = p.intermediateChordTo(q, 0.25)
            self.test('intermediateChordTo', i, '51.372294°N, 000.707192°E')
            self.test('intermediateChordTo', isinstance(i, LatLon), True)  # PYCHOK test attr?

            # courtesy of <https://GitHub.com/bakakaldsas>
            i = LatLon(52, 0, 100).intermediateChordTo(LatLon(48, 2, 200), 0.25)
            self.test('intermediateChordTo-h', i.height, '125.000', fmt='%.3f')  # PYCHOK test attr?

        if hasattr(LatLon, 'midpointTo'):
            m = p.midpointTo(q)
            self.test('midpointTo', m, '50.536327°N, 001.274614°E')  # PYCHOK test attr? # 50.5363°N, 001.2746°E

        if hasattr(LatLon, 'destination'):
            p = LatLon(51.4778, -0.0015)
            d = p.destination(7794, 300.7)
            self.test('destination', d, '51.513546°N, 000.098345°W' if Sph
                                   else '51.513526°N, 000.098038°W')  # 51.5135°N, 0.0983°W ???
            self.test('destination', d.toStr(F_DMS, 0), '51°30′49″N, 000°05′54″W' if Sph
                                                   else '51°30′49″N, 000°05′53″W')
            # <https://www.EdWilliams.org/avform.htm#LL>
            d = LAX.destination(100, 66, radius=R_NM) if Sph else LAX.destination(100, 66)
            self.test('destination', d.toStr(F_DM, prec=0), "34°37′N, 116°33′W" if Sph
                                                       else "33°57′N, 118°24′W")
            self.test('destination', d, '34.613647°N, 116.55116°W' if Sph
                                   else '33.950367°N, 118.399012°W')
            self.test('destination', d.toStr(F_RAD, prec=6), '0.604122N, 2.034201W' if Sph
                                                        else '0.592546N, 2.066453W')  # PYCHOK expected
            # <https://GeographicLib.SourceForge.io/html/python/examples.html>
            d = LatLon(-32.06, -115.74).destination(20e6, 225).toStr(F_D, prec=8)
            self.test('destination', d, '31.96383509°N, 064.37329146°E' if Sph
                                   else '32.11195529°N, 063.95925278°E', known=True)  # PYCHOK test attr?

        if hasattr(LatLon, 'alongTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.alongTrackDistanceTo(s, 96, TypeError.__name__, known=True)
                self.test('alongTrackDistanceTo', d, 62331.59, fmt='%.2f')  # 62331
            except TypeError as x:
                self.test('alongTrackDistanceTo', str(x), 'incompatible ...', known=True)  # PYCHOK test attr?
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, 62331.58, fmt='%.2f')

            # <https://www.EdWilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.alongTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('alongTrackDistanceTo', d, 99.588, fmt='%.3f')  # NM

            # courtesy of Rimvydas Naktinis
            p = LatLon(53.36366, -1.83883)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, -7702.7, fmt='%.1f')

            p = LatLon(53.35423, -1.60881)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, 7587.6, fmt='%.1f')  # PYCHOK test attr?

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.crossTrackDistanceTo(s, 96)
                self.test('crossTrackDistanceTo', d, TypeError.__name__, known=True)
            except TypeError as x:
                self.test('crossTrackDistanceTo', str(x), 'incompatible ...', known=True)  # PYCHOK test attr?
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, -307.55, fmt='%.2f')  # -307.5

            # <https://www.EdWilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.crossTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('crossTrackDistanceTo', d, 7.4524, fmt='%.4f')  # PYCHOK test attr? # XXX 7.4512 NM

        if hasattr(LatLon, 'greatCircle'):
            p = LatLon(53.3206, -1.7297)
            gc = p.greatCircle(96.0)
            self.test('greatCircle', gc, '(-0.79408, 0.12856, 0.59406)')  # PYCHOK test attr?

        p = LatLon(53.3206, -1.7297)
        q = LatLon(53.1887, 0.1334)
        self.test('cosineAndoyerLambertTo', p.cosineAndoyerLambertTo(q), '124801.098' if Sph else '125205.962', fmt='%.3f')
        self.test('cosineAndoyerLambertTo', q.cosineAndoyerLambertTo(p), '124801.098' if Sph else '125205.962', fmt='%.3f')

        self.test('cosineForsyheAndoyerLambertTo', p.cosineForsytheAndoyerLambertTo(q), '124801.098' if Sph else '125205.965', fmt='%.3f')
        self.test('cosineForsyheAndoyerLambertTo', q.cosineForsytheAndoyerLambertTo(p), '124801.098' if Sph else '125205.965', fmt='%.3f')

        self.test('cosineLawTo', p.cosineLawTo(q), '124801.098', fmt='%.3f')
        self.test('cosineLawTo', q.cosineLawTo(p), '124801.098', fmt='%.3f')

        self.test('equirectangularTo', p.equirectangularTo(q), '124804.754', fmt='%.3f')
        self.test('equirectangularTo', q.equirectangularTo(p), '124804.754', fmt='%.3f')

        self.test('euclideanTo', p.euclideanTo(q), '131273.287', fmt='%.3f')
        self.test('euclideanTo', q.euclideanTo(p), '131273.287', fmt='%.3f')

        self.test('flatLocalTo', p.flatLocalTo(q), '124804.754' if Sph else '125209.633', fmt='%.3f')
        self.test('flatLocalTo', q.flatLocalTo(p), '124804.754' if Sph else '125209.633', fmt='%.3f')

        self.test('flatPolarTo', p.flatPolarTo(q), '133663.257', fmt='%.3f')
        self.test('flatPolarTo', q.flatPolarTo(p), '133663.257', fmt='%.3f')

        self.test('haversineTo', p.haversineTo(q), '124801.098', fmt='%.3f')
        self.test('haversineTo', q.haversineTo(p), '124801.098', fmt='%.3f')

        self.test('hubenyTo', p.hubenyTo, p.flatLocalTo)
        self.test('hubenyTo', q.hubenyTo, q.flatLocalTo)

        self.test('thomasTo', p.thomasTo(q), '124801.098' if Sph else '125206.188', fmt='%.3f')
        self.test('thomasTo', q.thomasTo(p), '124801.098' if Sph else '125206.188', fmt='%.3f')

        self.test('vincentysTo', p.vincentysTo(q), '124801.098', fmt='%.3f')
        self.test('vincentysTo', q.vincentysTo(p), '124801.098', fmt='%.3f')

        if hasattr(LatLon, 'greatCircleTo'):
            gc = p.greatCircleTo(q)
            self.test('greatCircleTo', gc, '(-0.79408, 0.12859, 0.59406)')  # PYCHOK test attr?

        if isclockwise:
            f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            for _ in self.testiter():
                self.test('isclockwise', isclockwise(f), False)  # PYCHOK test attr?
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            for _ in self.testiter():
                self.test('isclockwise', isclockwise(t), True)  # PYCHOK test attr?
            for _ in self.testiter():
                try:
                    self.test('isclockwise', isclockwise(t[:2]), ValueError.__name__)
                except ValueError as x:
                    self.test('isclockwise', str(x), 'points (2): too few', known=True)  # PYCHOK test attr?
            # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
            p = LatLon(85, -135), LatLon(85, -45), LatLon(85, 45), LatLon(85, 135), LatLon(85, -135)
            for _ in self.testiter():
                try:
                    self.test('isclockwise', isclockwise(p), ValueError.__name__)  # PYCHOK test attr?
                except ValueError as x:
                    self.test('isclockwise', str(x), 'zero or polar area', known=True)  # PYCHOK test attr?

        if isconvex:
            f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
            for _ in self.testiter():
                self.test('isconvex', isconvex(f), False)  # PYCHOK test attr?
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            for _ in self.testiter():
                self.test('isconvex', isconvex(t), True)  # PYCHOK test attr?
            for _ in self.testiter():
                try:
                    self.test('isconvex', isconvex(t[:2]), ValueError.__name__)
                except ValueError as x:
                    self.test('isconvex', str(x), 'points (2): too few', known=True)  # PYCHOK test attr?

        if isenclosedBy:
            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            for _ in self.testiter():
                self.test('isenclosedBy1', isenclosedBy(LatLon(45.5, 1.5), b), True)  # PYCHOK test attr?
            for _ in self.testiter():  # on polygon point is outside
                self.test('isenclosedBy2', isenclosedBy((46, 2), b), False)  # PYCHOK test attr?
            for _ in self.testiter():  # on polygon point is outside
                self.test('isenclosedBy3', isenclosedBy((46, 1), b), False)  # PYCHOK test attr?
            for _ in self.testiter():  # on polygon edge is outside
                self.test('isenclosedBy4', isenclosedBy((46, 1.5), b), False)  # PYCHOK test attr?
            for _ in self.testiter():  # on polygon is False
                self.test('isenclosedBy5', isenclosedBy((45.5, 1), b), False)  # PYCHOK test attr?
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90), LatLon(85, -180)
            for _ in self.testiter():
                self.test('isenclosedBy6', isenclosedBy((90, 0), p), 'True')  # PYCHOK test attr?
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90)
            for _ in self.testiter():
                self.test('isenclosedBy7', isenclosedBy((90, 0), p), 'True')  # PYCHOK test attr?

        if hasattr(LatLon, 'initialBearingTo'):
            b = LHR.initialBearingTo(FRA)
            self.test('initialBearingTo', b, 102.432182 if Sph else 102.392291, fmt='%.6f')  # PYCHOK test attr?
        a = LHR.compassAngleTo(FRA, adjust=False)
        self.test('compassAngleTo', a, 100.017, fmt='%.3f')
        a = LHR.compassAngleTo(FRA)  # adjust=True
        self.test('compassAngleTo', a, 105.599, fmt='%.3f')

        if hasattr(LatLon, 'initialBearingTo'):
            b = FRA.initialBearingTo(LHR)
            self.test('initialBearingTo', b, 288.715918 if Sph else 288.676039, fmt='%.6f')  # PYCHOK test attr?
        a = FRA.compassAngleTo(LHR, adjust=False)
        self.test('compassAngleTo', a, 280.017, fmt='%.3f')
        a = FRA.compassAngleTo(LHR)  # adjust=True
        self.test('compassAngleTo', a, 285.599, fmt='%.3f')

        d = LHR.equirectangularTo(FRA)
        self.test('equirectangularTo', m2km(d), 592.185, fmt='%.3f')
        if hasattr(LatLon, 'distanceTo'):
            d = LHR.distanceTo(FRA)
            self.test('distanceTo', m2km(d), 591.831 if Sph else (591.831 if Nv else 593.571), fmt='%.3f')  # PYCHOK test attr?

        p = LatLon(0, 0)
        for a, b, d in ((1, 0, 0.0), ( 1, 1,  45.0), ( 0,  1,  90.0),
                                     (-1, 0, 180.0), (-1, -1, 225.0),
                                     (1, -1, 315.0), ( 0, -1, 270.0),
                                     (90, -1, 359.4)):
            q = LatLon(p.lat + a, p.lon + b)
            if not Nv:
                b = p.initialBearingTo(q)
                self.test('bearingTo', b, d, fmt='%.1f', known=True)
            c = p.compassAngleTo(q, adjust=False)
            self.test('compassAngleTo', c, d, fmt='%.1f')

        p.lat, p.lon = 52, '0'  # coverage
        q = LatLon(p.lat + 1, p.lon + 1)  # 45 unadjusted
        self.test('latlon2', q.latlon2(1), '(53.0, 1.0)')
        self.test('philam2', q.philam2(2), '(0.93, 0.02)')
        if not Nv:
            b = p.initialBearingTo(q)
            self.test('bearingTo', b, 31.0, fmt='%.0f')
        c = p.compassAngleTo(q, adjust=True)
        self.test('compassAngleTo', c, 31.0, fmt='%.0f')
        c = p.compassAngleTo(q, adjust=False)
        self.test('compassAngleTo', c, 45.0, fmt='%.0f')

        # check return types
        self.testReturnType(p.boundsOf(2, 4),       Bounds2Tuple, 'boundsOf')
        self.testReturnType(p.latlon,               LatLon2Tuple, 'latlon')
        self.testReturnType(p.latlon2(),            LatLon2Tuple, 'latlon2')
        self.testReturnType(p.latlonheight,         LatLon3Tuple, 'latlonheight')
        self.testReturnType(p.isequalTo(p),         bool,         'isequalTo')
        self.testReturnType(p.philam,               PhiLam2Tuple, 'philam')
        self.testReturnType(p.philamheight,         PhiLam3Tuple, 'philamheight')
        self.testReturnType(p.to2ab(),              PhiLam2Tuple, 'to2ab')
        self.testReturnType(p.to3llh(0),            LatLon3Tuple, 'to3llh')
        self.testReturnType(p.to3xyz(),             Vector3Tuple, 'to3xyz')
        self.testReturnType(p.xyz,                  Vector3Tuple, 'xyz')
        self.testReturnType(p.xyzh,                 Vector4Tuple, 'xyzh')

        self.testReturnType(p.compassAngleTo(q),    float, 'compassAngleTo')
        self.testReturnType(p.cosineLawTo(q),       float, 'cosineLawTo')
        self.testReturnType(p.euclideanTo(q),       float, 'euclideanTo')
        self.testReturnType(p.flatLocalTo(q),       float, 'flatLocalTo')
        self.testReturnType(p.flatPolarTo(q),       float, 'flatPolarTo')
        self.testReturnType(p.haversineTo(q),       float, 'haversineTo')
        self.testReturnType(p.hubenyTo(q),          float, 'hubenyTo')
        self.testReturnType(p.vincentysTo(q),       float, 'vincentysTo')

        if not Nv:
            # XXX prec=5 for NvectorBase.toStr vs prec=6 for Vector3Tuple.toStr
            self.test('toNvector', p.toNvector(), '(0.61566, 0.0, 0.78801)' if Sph
                                             else '(0.615661, 0.0, 0.788011)')  # PYCHOK test
        self.test('toVector',   p.toVector(), '(0.615661, 0.0, 0.788011)')
        self.test('toVector3d', p.toVector3d(), '(0.61566, 0.0, 0.78801)')

    def testReturnType(self, inst, clas, name):
        self.test(name, type(inst), clas)  # type(inst).__name__ == clas.__name__


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)

    if geographiclib:
        from pygeodesy import ellipsoidalKarney
        t.testLatLon(ellipsoidalKarney, Sph=False, Nv=False)

    t.testLatLon(ellipsoidalNvector, Sph=False, Nv=True)
    t.testLatLon(ellipsoidalVincenty, Sph=False, Nv=False)

    t.testLatLon(sphericalNvector, Sph=True, Nv=True)
    t.testLatLon(sphericalTrigonometry, Sph=True, Nv=False)

    t.results()
    t.exit()
