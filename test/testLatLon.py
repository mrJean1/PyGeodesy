
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '21.09.30'

from base import GeodSolve, geographiclib, isPyPy, isPython2, TestsBase

from pygeodesy import F_D, F_DM, F_DMS, F_RAD, R_M, R_NM, \
                      degrees, fstr, Height, isclockwise, isconvex, \
                      isenclosedBy, isnear0, ispolar, m2km, m2NM, \
                      IntersectionError, VincentyError  # PYCHOK expected
from pygeodesy.namedTuples import Bounds2Tuple, \
                                  LatLon2Tuple, LatLon3Tuple, \
                                  PhiLam2Tuple, PhiLam3Tuple, \
                                  Vector3Tuple, Vector4Tuple


class Tests(TestsBase):

    def testLatLon(self, module, Sph=False, Nv=False, X=False, GS=False):  # MCCABE 108

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
        q = LatLon(*map(degrees, p.philam))
        self.test('isequalTo', q.isequalTo(p), True)

        self.test('latlon2', fstr(q.latlon2(5)), '52.20472, 0.14056')
        self.test('latlon2', fstr(q.latlon2(4)), '52.2047, 0.1406')
        self.test('latlon2', fstr(q.latlon2(3)), '52.205, 0.141')
        self.test('latlon2', fstr(q.latlon2(2)), '52.2, 0.14')
        self.test('latlon2', fstr(q.latlon2(1)), '52.2, 0.1')
        self.test('latlon2', fstr(q.latlon2()),  '52.0, 0.0')

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

        # <https://www.Numericana.com/answer/formula.htm#geodetic>
        p = LatLon( 42.96125,  -85.655719, height=195)  # Grand Rapids, MI
        q = LatLon(-37.813611, 144.963056, height=31)  # Melbourne, AU
        self.test('chordTo', q.chordTo(p),           '12036677.26' if Sph else '12029263.15', fmt='%.2f')
        self.test('chordTo', p.chordTo(q, height=0), '12036463.78' if Sph else '12029049.69', fmt='%.2f')

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
            self.test('initialBearingTo', b, '156.1666' if Sph else '156.1106', fmt='%.4f')  # 156.2
            b = LAX.initialBearingTo(JFK)
            self.test('initialBearingTo', b, '65.8921' if Sph else '65.9335', fmt='%.4f')
            # courtesy Inmars Didaitis (bakakaldsas) <https://GitHub.com/mrJean1/PyGeodesy/issues/56>
            a, b = LatLon(5, 5), LatLon(10, 5)  # same lon
            self.test('initialBearingTo', a.initialBearingTo(b),   0.0, prec=1)
            self.test('initialBearingTo', b.initialBearingTo(a), 180.0, prec=1)

        if hasattr(LatLon, 'finalBearingTo'):
            b = p.finalBearingTo(q)
            self.test('finalBearingTo', b, '157.8904' if Sph else '157.8345', fmt='%.4f')
            b = LAX.finalBearingTo(JFK)
            self.test('finalBearingTo', b, '93.8581' if Sph else '93.9034', fmt='%.4f')

        if hasattr(LatLon, 'bearingTo2'):
            b = p.bearingTo2(q)
            self.test('bearingTo2', fstr(b, prec=4), '156.1666, 157.8904' if Sph else '156.1106, 157.8345')  # 156.2
            # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
            b = LatLon(85, -135), LatLon(85, -45), LatLon(85, 45), LatLon(85, 135), LatLon(85, -135)
            self.test('ispolar', ispolar(b), True)

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
            self.test('distanceTo', d, 2145 if Sph else (3972863 if Nv else 3981601), fmt='%.0f')
            if not Nv:  # <https://GeographicLib.SourceForge.io/html/python/examples.html>
                self.test('antipodal', WNZ.isantipodeTo(SAL, eps=0.1), False)
                try:
                    d = WNZ.distanceTo(SAL, wrap=False)
                    self.test('distanceTo dateline', d, '19119590.551' if Sph else '19959679.267', fmt='%.3f', known=True)
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo dateline', x, 'no convergence ...', known=True)
                try:
                    d = WNZ.distanceTo(SAL, wrap=True)
                    self.test('distanceTo unrolled', d, '19119590.551' if Sph else '19959679.267', fmt='%.3f', known=True)
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo unrolled', x, 'no convergence ...', known=True)
                self.test('antipodal', BJS.isantipodeTo(SFO, eps=0.1), False)
                d = BJS.distanceTo(SFO, wrap=False)
                self.test('distanceTo dateline', d, '9491735' if Sph else '9513998', fmt='%.0f')
                d = BJS.distanceTo(SFO, wrap=True)
                self.test('distanceTo unrolled', d, '9491735' if Sph else '9513998', fmt='%.0f')

            # <https://GitHub.com/ChrisVeness/geodesy/issues/64>
            d = LatLon(20, 0).distanceTo(LatLon(-2, 180))
            self.test('distanceTo', d, '18013602.92' if Sph or Nv else ('18003740.39' if X else '18012714.66'), fmt='%.2f')
            try:
                d = LatLon(0, 0).distanceTo(LatLon(0, 180))  # antipodal
                self.test('distanceTo', d, '20015114.35' if Sph else '20003931.46', fmt='%.2f', known=Nv or X or GS)  # PYCHOK 0.0 for Nv ...
            except ValueError as x:
                x = str(x)
                self.test('distanceTo', x, 'ambiguous, antipodal ...', known=True)

        if hasattr(LatLon, 'distanceTo3'):  # and not Sph:
            for w in (False, True):
                try:
                    d = WNZ.distanceTo3(SAL, wrap=w)  # XXX expected values?
                    self.test('distanceTo3 dateline', fstr(d, prec=4), '19959679.2674, 161.0677, 18.8252', known=True)
                except VincentyError as x:
                    x = str(x)
                    self.test('distanceTo3 dateline', x, 'no convergence ...', known=True)
                d = BJS.distanceTo3(SFO, wrap=w)
                self.test('distanceTo3 dateline', fstr(d, prec=4), '9513997.9901, 42.9164, 138.8903')

        if hasattr(LatLon, 'intermediateTo'):
            i = p.intermediateTo(q, 0.25)
            self.test('intermediateTo', i, '51.372084°N, 000.707337°E' if Sph
                                     else ('51.372294°N, 000.707192°E' if Nv
                                     else  '51.372275°N, 000.707253°E'))
            self.test('intermediateTo', isinstance(i, LatLon), True)

            if hasattr(p, 'distanceTo'):
                d = p.distanceTo(q)
                self.test('intermediateTo', d, '404279.721' if Sph or Nv
                                          else '404607.806', prec=3)

            i = p.intermediateTo(q, 5)
            self.test('intermediateTo+5', i, '35.160975°N, 008.989542°E' if Sph
                                       else ('35.560239°N, 008.833512°E' if Nv
                                       else  '35.139582°N, 008.994368°E'))
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo+5', p.distanceTo(i) / d, '4.885' if Nv and not Sph else '5.000', prec=3)

            i = p.intermediateTo(q, -4)
            self.test('intermediateTo-4', i, '64.911647°N, 013.726301°W' if Sph
                                       else ('64.570387°N, 013.156352°W' if Nv
                                       else  '64.894124°N, 013.705689°W'))
            if hasattr(p, 'distanceTo'):
                self.test('intermediateTo-4', p.distanceTo(i) / d, '3.885' if Nv and not Sph else '4.000', prec=3)

            # courtesy of <https://GitHub.com/bakakaldsas>
            i = LatLon(52, 0, 100).intermediateTo(LatLon(48, 2, 200), 0.25)
            self.test('intermediateTo-h', i.height, '125.000', prec=3)

        if hasattr(LatLon, 'intermediateChordTo'):
            i = p.intermediateChordTo(q, 0.25)
            self.test('intermediateChordTo', i, '51.372294°N, 000.707192°E')
            self.test('intermediateChordTo', isinstance(i, LatLon), True)

            # courtesy of <https://GitHub.com/bakakaldsas>
            i = LatLon(52, 0, 100).intermediateChordTo(LatLon(48, 2, 200), 0.25)
            self.test('intermediateChordTo-h', i.height, '125.000', fmt='%.3f')

        if hasattr(LatLon, 'midpointTo'):
            m = p.midpointTo(q)
            self.test('midpointTo', m, '50.536327°N, 001.274614°E')  # 50.5363°N, 001.2746°E

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
                                   else '32.11195529°N, 063.95925278°E', known=True)

        if hasattr(LatLon, 'alongTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.alongTrackDistanceTo(s, 96, TypeError.__name__)
                self.test('alongTrackDistanceTo', d, '62331.59', fmt='%.2f', known=True)  # 62331
            except TypeError as x:
                self.test('alongTrackDistanceTo', str(x), 'incompatible ...', known=True)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '62331.58', fmt='%.2f')

            # <https://www.EdWilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.alongTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('alongTrackDistanceTo', d, '99.588', fmt='%.3f')  # NM

            # courtesy of Rimvydas Naktinis
            p = LatLon(53.36366, -1.83883)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '-7702.7', fmt='%.1f')

            p = LatLon(53.35423, -1.60881)
            d = p.alongTrackDistanceTo(s, e)
            self.test('alongTrackDistanceTo', d, '7587.6', fmt='%.1f')

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            s = LatLon(53.3206, -1.7297)
            e = LatLon(53.1887, 0.1334)
            p = LatLon(53.2611, -0.7972)
            try:
                d = p.crossTrackDistanceTo(s, 96)
                self.test('crossTrackDistanceTo', d, TypeError.__name__, known=True)
            except TypeError as x:
                self.test('crossTrackDistanceTo', str(x), 'incompatible ...', known=True)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', fmt='%.2f')  # -307.5

            # <https://www.EdWilliams.org/avform.htm#XTE>
            p = LatLon(34.5, -116.5)  # 34:30N, 116:30W
            d = p.crossTrackDistanceTo(LAX, JFK, radius=Rav)
            self.test('crossTrackDistanceTo', d, '7.4524', fmt='%.4f')  # XXX 7.4512 NM

        if hasattr(LatLon, 'greatCircle'):
            p = LatLon(53.3206, -1.7297)
            gc = p.greatCircle(96.0)
            self.test('greatCircle', gc, '(-0.79408, 0.12856, 0.59406)')

        if hasattr(LatLon, 'nearestOn6'):
            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            p = LatLon(1, 1).nearestOn6(b, height=0)
            self.test('neareston6', p, '(LatLon(45°00′00.0″N, 001°00′00.0″E), 4755443.4294, 0.0, 1, LatLon(45°00′00.0″N, 001°00′00.0″E), LatLon(45°00′00.0″N, 001°00′00.0″E))', known=not X)
            p = LatLon(45.5, 2.5).nearestOn6(b, height=0)
            self.test('neareston6', p, '(LatLon(45°30′03.94″N, 002°00′00.0″E), 39078.729285, 1.501072, 2, LatLon(45°00′00.0″N, 002°00′00.0″E), LatLon(46°00′00.0″N, 002°00′00.0″E))', known=not X)

        p = LatLon(53.3206, -1.7297)
        q = LatLon(53.1887, 0.1334)
        self.test('chordTo', p.chordTo(q), '124799.103' if Sph else '125203.963', fmt='%.3f')

        self.test('cosineAndoyerLambertTo', p.cosineAndoyerLambertTo(q), '124801.098' if Sph else '125205.962', fmt='%.3f')
        self.test('cosineAndoyerLambertTo', q.cosineAndoyerLambertTo(p), '124801.098' if Sph else '125205.962', fmt='%.3f')

        self.test('cosineForsyheAndoyerLambertTo', p.cosineForsytheAndoyerLambertTo(q), '124801.098' if Sph else '125205.965', fmt='%.3f')
        self.test('cosineForsyheAndoyerLambertTo', q.cosineForsytheAndoyerLambertTo(p), '124801.098' if Sph else '125205.965', fmt='%.3f')

        self.test('cosineLawTo', p.cosineLawTo(q), '124801.098', fmt='%.3f')
        self.test('cosineLawTo', q.cosineLawTo(p), '124801.098', fmt='%.3f')

        self.test('equirectangularTo', p.equirectangularTo(q), '124804.754', fmt='%.3f')
        self.test('equirectangularTo', q.equirectangularTo(p), '124804.754', fmt='%.3f')

        self.test('euclideanTo', p.euclideanTo(q), '130015.089', prec=3)  # XXX 131273.287 @ _0_5
        self.test('euclideanTo', q.euclideanTo(p), '130015.089', prec=3)

        self.test('flatLocalTo', p.flatLocalTo(q), '124804.754' if Sph else '125209.633', fmt='%.3f')
        self.test('flatLocalTo', q.flatLocalTo(p), '124804.754' if Sph else '125209.633', fmt='%.3f')

        self.test('flatPolarTo', p.flatPolarTo(q), '133663.257', fmt='%.3f')
        self.test('flatPolarTo', q.flatPolarTo(p), '133663.257', fmt='%.3f')

        self.test('hartzell', p.hartzell(), p)  # '53.3206°N, 001.7297°W'
        c = p.toCartesian()  # pos 1,000 km up, los non-None, height is near 0
        t = c.toLatLon(height=1e6).hartzell(los=c.negate()).toStr(form=F_D, prec=6, m=None)
        self.test('hartzell', t, '53.3206°N, 001.7297°W' if Sph else '53.349541°N, 001.7297°W')

        h = p.height4().h
        self.test('height4', h, float(p.height), known=abs(h) < 2e-9)
        self.test('height4', p.height4(earth=R_M).toStr(prec=1), '(3803904.2, -114870.8, 5109488.3, 0.0)' if Sph
                                                            else '(3820333.9, -115367.0, 5097204.4, -6584.9)')
        self.test('height4', p.height4(LatLon=LatLon, height=0).toStr(prec=1), '53°19′14.2″N, 001°43′46.9″W')

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
            self.test('greatCircleTo', gc, '(-0.79408, 0.12859, 0.59406)')

        if isclockwise:
            f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            for _ in self.testiter():
                self.test('isclockwise', isclockwise(f), False)
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            for _ in self.testiter():
                self.test('isclockwise', isclockwise(t), True)
            for _ in self.testiter():
                try:
                    self.test('isclockwise', isclockwise(t[:2]), ValueError.__name__)
                except ValueError as x:
                    self.test('isclockwise', str(x), 'points (2): too few', known=True)
            # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
            p = LatLon(85, -135), LatLon(85, -45), LatLon(85, 45), LatLon(85, 135), LatLon(85, -135)
            for _ in self.testiter():
                try:
                    self.test('isclockwise', isclockwise(p), ValueError.__name__)
                except ValueError as x:
                    self.test('isclockwise', str(x), 'zero or polar area', known=True)

        if isconvex:
            f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
            for _ in self.testiter():
                self.test('isconvex', isconvex(f), False)
            t = LatLon(45,1), LatLon(46,1), LatLon(46,2), LatLon(45,1)
            for _ in self.testiter():
                self.test('isconvex', isconvex(t), True)
            for _ in self.testiter():
                try:
                    self.test('isconvex', isconvex(t[:2]), ValueError.__name__)
                except ValueError as x:
                    self.test('isconvex', str(x), 'points (2): too few', known=True)

        if isenclosedBy:
            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            for _ in self.testiter():
                self.test('isenclosedBy1', isenclosedBy(LatLon(45.5, 1.5), b), True)
            for _ in self.testiter():  # on polygon point is outside
                self.test('isenclosedBy2', isenclosedBy((46, 2), b), False)
            for _ in self.testiter():  # on polygon point is outside
                self.test('isenclosedBy3', isenclosedBy((46, 1), b), False)
            for _ in self.testiter():  # on polygon edge is outside
                self.test('isenclosedBy4', isenclosedBy((46, 1.5), b), False)
            for _ in self.testiter():  # on polygon is False
                self.test('isenclosedBy5', isenclosedBy((45.5, 1), b), False)
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90), LatLon(85, -180)
            for _ in self.testiter():
                self.test('isenclosedBy6', isenclosedBy((90, 0), p), 'True')
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90)
            for _ in self.testiter():
                self.test('isenclosedBy7', isenclosedBy((90, 0), p), 'True')

        if hasattr(LatLon, 'isenclosedBy'):
            # courtesy Maranov <https://GitHub.com/mrJean1/PyGeodesy/issues/53>
            b = LatLon(0, 0), LatLon(0, 1), LatLon(0, 1), LatLon(1, 1), LatLon(1, 0)
            for t in ('CCW', 'CW ', 'CCW'):
                self.test('isenclosedBy-' + t, LatLon(0.5, 0.5).isenclosedBy(b), True)
                b = tuple(reversed(b))

        if hasattr(LatLon, 'initialBearingTo'):
            b = LHR.initialBearingTo(FRA)
            self.test('initialBearingTo', b, '102.432182' if Sph else '102.392291', fmt='%.6f')
        a = LHR.compassAngleTo(FRA, adjust=False)
        self.test('compassAngleTo', a, '100.017', fmt='%.3f')
        a = LHR.compassAngleTo(FRA)  # adjust=True
        self.test('compassAngleTo', a, '105.599', fmt='%.3f')

        if hasattr(LatLon, 'initialBearingTo'):
            b = FRA.initialBearingTo(LHR)
            self.test('initialBearingTo', b, '288.715918' if Sph else '288.676039', fmt='%.6f')
        a = FRA.compassAngleTo(LHR, adjust=False)
        self.test('compassAngleTo', a, '280.017', fmt='%.3f')
        a = FRA.compassAngleTo(LHR)  # adjust=True
        self.test('compassAngleTo', a, '285.599', fmt='%.3f')

        d = LHR.equirectangularTo(FRA)
        self.test('equirectangularTo', m2km(d), '592.185', fmt='%.3f')
        if hasattr(LatLon, 'distanceTo'):
            d = LHR.distanceTo(FRA)
            self.test('distanceTo', m2km(d), '591.831' if Sph else ('591.831' if Nv else '593.571'), fmt='%.3f')

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
#       self.testReturnType(p.height4(),            Vector4Tuple, 'height4')
        self.testReturnType(p.isequalTo(p),         bool,         'isequalTo')
        self.testReturnType(p.philam,               PhiLam2Tuple, 'philam')
        self.testReturnType(p.philamheight,         PhiLam3Tuple, 'philamheight')
#       self.testReturnType(p.to2ab(),              PhiLam2Tuple, 'to2ab')
#       self.testReturnType(p.to3llh(0),            LatLon3Tuple, 'to3llh')
#       self.testReturnType(p.to3xyz(),             Vector3Tuple, 'to3xyz')
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
                                             else '(0.615661, 0.0, 0.788011)')
        self.test('toVector',   p.toVector(), '(0.615661, 0.0, 0.788011)')
        self.test('toVector3d', p.toVector3d(), '(0.61566, 0.0, 0.78801)')

        # courtesy AleixDev <https://GitHub.com/mrJean1/PyGeodesy/issues/43>
        def _known(p):
            return int(p.lat) == 42 and int(p.lon) == 2

        n  = 'trilaterate5 (%s) .' % (LatLon.__module__,)
        d  = 5110  # meter
        p1 = LatLon(42.688839, 2.438857)
        p2 = LatLon(42.635421, 2.522570)

        p3 = LatLon(42.630788, 2.500713)
        if Nv:  # coverage
            t = p1.trilaterate5(d, p2, d, p3, d, area=False, eps=250)  # intersection
            self.test(n + 'min', t.min, '223.305', prec=3, known=120 < t.min < 150)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.67456065°N, 002.49539502°E', known=_known(p))
            self.test(n + 'min- is .maxPoint', t.minPoint is t.maxPoint, True, known=Sph)
            self.test(n + 'n', t.n, t.n)

        try:
            t = p1.trilaterate5(d, p2, d, p3, d, area=True)  # overlap
            self.test(n + 'min', t.min, '313.671' if Sph
                                  else ('305.091' if geographiclib
                                  else  '311.234'), prec=3, known= 300 < t.min <  320)
            p = t.minPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.66937229°N, 002.48639477°E' if Sph
                                                   else ('42.66933643°N, 002.48620262°E' if geographiclib
                                                   else  '42.66938776°N, 002.48641176°E'), known=_known(p))
            self.test(n + 'max', t.max, '1591.044' if Sph
                                  else ('1592.545' if geographiclib
                                  else  '1586.951'), prec=3, known=1570 < t.max < 1610)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.65153054°N, 002.46822157°E' if Sph
                                                   else ('42.65141232°N, 002.46816989°E' if geographiclib
                                                   else  '42.65153156°N, 002.46821899°E'), known=_known(p))
            self.test(n + 'n', t.n, t.n)

            t = p1.trilaterate5(d, p2, d, p3, d, area=False, eps=250)  # intersection
            self.test(n + 'min', t.min, '133.815' if Sph
                                  else ('127.229' if geographiclib
                                  else  '137.897'), prec=3, known=120 < t.min < 150)
            self.test(n + 'max', t.max, '160.242' if Sph
                                  else ('152.612' if geographiclib
                                  else  '148.175'), prec=3, known=120 < t.max < 170)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.67817811°N, 002.49966641°E' if Sph
                                                   else ('42.67815375°N, 002.49950041°E' if geographiclib
                                                   else  '42.67811504°N, 002.49959193°E'), known=_known(p))
            self.test(n + 'n', t.n, t.n)

            p3 = LatLon(42.64540, 2.504811)
            t = p1.trilaterate5(d, p2, d, p3, d, area=True)  # overlap
            x = '2403.293' if Sph else ('2400.293' if geographiclib else '2399.908')
            self.test(n + 'min', t.min, x, prec=3, known=2380 < t.min < 2420)
            self.test(n + 'max', t.max, x, prec=3, known=2380 < t.max < 2420)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.66135649°N, 002.47981645°E' if Sph
                                                   else ('42.66128984°N, 002.47973818°E' if geographiclib
                                                   else  '42.6613586°N, 002.47981223°E'), known=_known(p))
            self.test(n + 'min- is .maxPoint', t.minPoint is t.maxPoint, True)
            self.test(n + 'n', t.n, t.n)

            t = p1.trilaterate5(d, p2, d, p3, d, area=False, eps=1500)  # intersection
            self.test(n + 'min', t.min, '1340.608' if Sph
                                  else ('1343.743' if geographiclib
                                  else  '1332.749'), prec=3, known=1320 < t.min < 1360)
            p = t.minPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.69128229°N, 002.50129001°E' if Sph
                                                   else ('42.69131964°N, 002.50112167°E' if geographiclib
                                                   else  '42.69124153°N, 002.50124031°E'), known=_known(p))
            self.test(n + 'max', t.max, '1499.220' if Sph
                                  else ('1445.554' if geographiclib
                                  else  '1450.709'), prec=3, known=1420 < t.max < 1520)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.64295864°N, 002.44242391°E' if Sph
                                                   else ('42.67815375°N, 002.49950041°E' if geographiclib
                                                   else  '42.67811504°N, 002.49959193°E'), known=_known(p))
            self.test(n + 'n', t.n, t.n)

            k = (isPyPy and isPython2) or Sph or not X
            t = repr(p1.radii11(p2, p3))
            self.test('radii11', t, t, known=k)

            n = 'circum3 (%s) .' % (LatLon.__module__,)
            t = p1.circum3(p2, p3)
            c = t.center
            self.test(n + 'radius', t.radius, '57792.067', prec=3, known=k)
            self.test(n + 'center', c, '43.053532°N, 002.943255°E, -261.66m', known=k)
            d = t.deltas
            self.test(n + 'deltas', d.toStr(prec=3), '(0.0, 0.0, 11.383)', known=k or isnear0(d.lat, 1e-8) or isnear0(d.lon, 1e-8))

            c.height = 0
            self.test(n + 'd1', p1.distanceTo(c), '57792.858', prec=3, known=k)
            self.test(n + 'd2', p2.distanceTo(c), '57792.859', prec=3, known=k)
            self.test(n + 'd3', p3.distanceTo(c), '57792.859', prec=3, known=k)

            self.test(n + 'datum', c.datum, p1.datum)
            self.test(n + 'Ecef',  c.Ecef,  p1.Ecef)

        except ImportError as x:
            self.skip(str(x), n=27)
        except IntersectionError as x:
            self.test(n + 'error', str(x), str(x))
        except NotImplementedError as x:
            self.test(n + 'error', str(x), str(x) if Nv else NotImplementedError.__name__)

        try:  # need numpy
            k = (isPyPy and isPython2) or Sph or not X
            t = repr(p1.radii11(p2, p3))
            self.test('radii11', t, t, known=k)

            n = 'circum4 (%s) .' % (LatLon.__module__,)
            t = p1.circum4_(p2, p3)
            c = t.center
            self.test(n + 'radius', t.radius, '3184256.748', prec=3, known=k)
            self.test(n + 'center', c, '43.054367°N, 002.942573°E, -3183993.92m', known=k)
            self.test(n + 'rank', t.rank, 3, known=k)
            self.test(n + 'residuals', t.residuals, (), known=k)

            c.height = 0
            self.test(n + 'd1', p1.distanceTo(c), '57818.033', prec=3, known=k)
            self.test(n + 'd2', p2.distanceTo(c), '57834.176', prec=3, known=k)
            self.test(n + 'd3', p3.distanceTo(c), '57830.992', prec=3, known=k)

            self.test(n + 'datum', c.datum, p1.datum)
            self.test(n + 'Ecef',  c.Ecef,  p1.Ecef)

            n = 'circin6 (%s) .' % (LatLon.__module__,)
            t = LatLon(0, 0).radii11(LatLon(10, 0), LatLon(0, 10))
            r = t.toRepr()
            self.test('radii11', r, r, known=k)
            self.test(n + 'rB+rC', t.rB + t.rC, t.a, known=k, prec=3)
            self.test(n + 'rC+pA', t.rC + t.rA, t.b, known=k, prec=3)
            self.test(n + 'rA+rB', t.rA + t.rB, t.c, known=k, prec=3)

            t = LatLon(0, 0).circin6(LatLon(10, 0), LatLon(0, 10))  # XXX no intersection
            c = t.center
            self.test(n + 'radius', t.radius, '325058.721', prec=3, known=k)
            self.test(n + 'center', c, '02.948531°N, 002.932537°E, -40041.19m', known=k)
            self.test(n + 'deltas', t.deltas, '(0.0, 0.0, 0.090491)', known=True)  # XXX

            self.test(n + 'cA', t.cA, '05.04314°N, 005.014578°E, -48104.09m', known=k)
            self.test(n + 'cB', t.cB, '00.0°N, 002.941713°E, -20168.62m', known=k)
            self.test(n + 'cC', t.cC, '02.961566°N, 000.0°E, -20113.46m', known=k)

            # c.height = 0
            self.test(n + 'dA', c.distanceTo(t.cA), '327263.596', prec=3, known=k)
            self.test(n + 'dB', c.distanceTo(t.cB), '326036.153', prec=3, known=k)
            self.test(n + 'dC', c.distanceTo(t.cC), '326020.432', prec=3, known=k)

        except ImportError as x:
            self.skip(str(x), n=22)

    def testReturnType(self, inst, clas, name):
        self.test(name, type(inst), clas)  # type(inst).__name__ == clas.__name__


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)

    t.testLatLon(sphericalNvector, Sph=True, Nv=True)
    t.testLatLon(sphericalTrigonometry, Sph=True)

    t.testLatLon(ellipsoidalNvector, Nv=True)
    t.testLatLon(ellipsoidalVincenty)

    if geographiclib:
        from pygeodesy import ellipsoidalKarney
        t.testLatLon(ellipsoidalKarney)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testLatLon(ellipsoidalGeodSolve, GS=True)

    t.testLatLon(ellipsoidalExact, X=True)

    t.results()
    t.exit()
