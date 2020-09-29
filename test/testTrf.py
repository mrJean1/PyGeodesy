
# -*- coding: utf-8 -*-

u'''Test Terrestrial Reference Frame (TRF) implementation.

All tests transcribed from Chris Veness (C) 2006-2019 U{Geodesy tools for conversions between
reference frames<https://www.Movable-Type.co.UK/scripts/geodesy-library.html>} JavaScript original
<https://GitHub.com/chrisveness/geodesy/blob/master/test/latlon-ellipsoidal-referenceframe-tests.js>
'''

__all__ = ('Tests',)
__version__ = '20.09.27'

from base import TestsBase

from pygeodesy import date2epoch, Epoch, epoch2date, F_D, F_DMS, RefFrames, TRFError


class Tests(TestsBase):

    def testTrf(self, Cartesian, LatLon, ellipsoidal):  # MCCABE 17

        self.test('ellipsoidal' + ellipsoidal, '...', '...')

        p = LatLon(51.47788, -0.00147, reframe=RefFrames.ITRF2000)
        x = p.convertRefFrame(RefFrames.ETRF2000)
        self.test('convertRefFrame', x.toStr(form=F_D, prec=8), '51.47787826°N, 000.00147125°W', known=True)

        c = Cartesian(4027893.924, 307041.993, 4919474.294)
        x = c.toLatLon()
        x.reframe = RefFrames.ITRF2000
        self.testCopy(x.reframe)
        self.test('toLatLon', x.toStr(form=F_D, prec=4), '50.7978°N, 004.3592°E, +148.96m')
        c = Cartesian(3980574.247, -102.127, 4966830.065)
        x = c.convertRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2000)
        self.test('convertRefFrame', x, '[3980574.395, -102.214, 4966829.941]')

        p = LatLon(0, 0, reframe=RefFrames.ITRF2000)
        x = p.convertRefFrame(RefFrames.ITRF2000)
        self.test('Nil', x is p, True)
        self.testCopy(x.reframe)
        c = Cartesian(1, 2, 3)
        x = c.convertRefFrame(RefFrames.ITRF2000, RefFrames.ITRF2000)
        self.test('Nil', x is c, True)

        p = LatLon(0, 0, reframe=RefFrames.NAD83)
        x = p.convertRefFrame(RefFrames.ITRF2014)  # # via ITRF2000
        self.test('reframe', x.reframe == RefFrames.ITRF2014, True)
        x = p.convertRefFrame(RefFrames.NAD83)
        self.test('Roundtrip', x == p, True)
        self.test('reframe', x.reframe == RefFrames.NAD83, True)
        self.testCopy(x.reframe)

        # Dawson, J. & Woods, A., Appendix A, Journal of Applied Geodesy 4 (2010)
        p = LatLon('23°40′12.41482″S', '133°53′7.86712″E', height=603.2562, reframe=RefFrames.ITRF2005, epoch=2010.4559)
        x = p.convertRefFrame(RefFrames.GDA94)
        self.test('Geodetic', x.toStr(form=F_DMS, prec=5), '23°40′12.44582″S, 133°53′07.84795″E, +603.34m')  # +603.3361m
        c = x.toCartesian()
        self.test('Cartesian', c.toStr(prec=4), '[-4052051.7614, 4212836.1945, -2545106.0147]')

        x = p.convertRefFrame(RefFrames.GDA94)  # epoch 2010.4559
        self.test('Geodetic', x.toStr(form=F_DMS, prec=5), '23°40′12.44582″S, 133°53′07.84795″E, +603.34m')  # +603.3361m'
        c = x.toCartesian()
        self.test('Cartesian', c.toStr(prec=4), '[-4052051.7614, 4212836.1945, -2545106.0147]')

        x = x.convertRefFrame(RefFrames.ITRF2005)  # epoch 2010.4559
        self.test('Roundtrip', x.toStr(form=F_DMS, prec=5), '23°40′12.41482″S, 133°53′07.86712″E, +603.26m')  # +603.2562m

        # <https://GitHub.com/OSGeo/proj.4/blob/2aaf53/test/gie/more_builtins.gie#L357>
        c = Cartesian(3370658.37800, 711877.31400, 5349787.08600)  # Proj4 Onsala observatory
        x = c.convertRefFrame(RefFrames.ITRF93, RefFrames.ITRF2000, 2017)
        self.test('GNSStrans', x.toStr(prec=5), '[3370658.18892, 711877.42369, 5349787.1243]')  # accurate to within 0.02 mm

        # <https://www.NGS.NOAA.gov/cgi-bin/ds_mark.prl?PidBox=kg0640>
        p = LatLon('39 13 26.71220', '098 32 31.74540', height=573.961, reframe=RefFrames.NAD83, epoch=2010.0)
        c = p.toCartesian()  # NGS Data Sheet Meades Ranch
        self.test('Cartesian', c, '[-734972.563, 4893188.492, 4011982.811]')

        # <https://EPNCB.OMA.BE/_productsservices/coord_trans> (tutorial)
        c = Cartesian(4027894.006, 307045.600, 4919474.910)
        x = c.convertRefFrame(RefFrames.ITRF91, RefFrames.ITRF2005, 2007)
        self.test('EUREF C1', x.toStr(prec=4), '[4027894.0444, 307045.6209, 4919474.8613]')
        x = c.convertRefFrame(RefFrames.ITRF91, RefFrames.ITRF2005, 2007)
        self.test('EUREF C2', x.toStr(prec=4), '[4027894.0444, 307045.6209, 4919474.8613]')
        x = c.convertRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2000, 2012)
        self.test('EUREF C4', x.toStr(prec=4), '[4027894.3559, 307045.2508, 4919474.6447]')
        x = c.convertRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2014, 2012)
        self.test('EUREF C5', x.toStr(prec=4), '[4027894.3662, 307045.253, 4919474.6263]')

        try:
            t = LatLon(0, 0, reframe='ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe) ('ITRF2000'): not a RefFrame")

        try:
            t = LatLon(0, 0, reframe=RefFrames.ITRF2000, epoch=1899)
        except TRFError as x:
            t = str(x)
        self.test(TRFError.__name__, t, "epoch (1899): below 1900 limit")

        try:
            t = LatLon(0, 0, reframe=RefFrames.ITRF2000).convertRefFrame('ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe2) ('ITRF2000'): not a RefFrame")

        try:
            t = LatLon(0, 0).convertRefFrame(RefFrames.ITRF2000)
        except TRFError as x:
            t = str(x)
        self.test('TRFError', t, 'no conversion: LatLon(00°00′00.0″N, 000°00′00.0″E).reframe MISSING')

        c = Cartesian(0, 0, 0)
        try:
            t = c.convertRefFrame('ITRF2000', RefFrames.ITRF2000)
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe2) ('ITRF2000'): not a RefFrame")

        try:
            t = c.convertRefFrame(RefFrames.ITRF2000, 'ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe) ('ITRF2000'): not a RefFrame")

#       try:
#           t = c.convertRefFrame(RefFrames.ITRF2000, RefFrames.ITRF2000, '2000')
#       except TypeError as x:
#           t = str(x)
#       self.test('TypeError', t, "type(epoch) ('2000'): not scalar")

    def testEpoch(self):

        try:  # coverage
            e = date2epoch(None, 1, 2)
            self.test('epoch', e, TRFError.__name__)
        except TRFError as x:
            t = str(x)
            self.test('TRFError', t, t)

        r = RefFrames.GDA94
        t = r.toStr()
        self.test('toStr', t, "name='GDA94', epoch=1994, ellipsoid=Ellipsoid(name='GRS80')")
        self.test('str', str(r),t)
        t = r.toStr2()
        self.test('toStr2', t, "RefFrame(name='GDA94', epoch=1994, ellipsoid=Ellipsoid(name='GRS80')")
        self.test('repr', repr(r), t)

        for y, m, d, x in ((2020,  1,  1, 2020.003),
                           (2020,  4,  1, 2020.251),
                           (2020,  7,  1, 2020.5),
                           (2020, 10,  1, 2020.751),
                           (2020, 12, 31, 2021.0)):
            e = date2epoch(y, m, d)
            self.test('epoch', e, x, fmt='%.3f')
            t = epoch2date(e)
            self.test('y-m-d', t, (y, m, d), known=t[0] == 2021)

        e = Epoch(Epoch=2020.)
        self.test(e.toRepr() + '.std_repr', e.std_repr, False)
        d = m = 0
        for n in (0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31):
            m += 1
            d += n
            e = 2020.001 + d / 366.0
            self.test(Epoch(Epoch=e).toRepr(), epoch2date(e), (2020, m, 1), known=m == 13)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalKarney as K, \
                          ellipsoidalNvector as N, \
                          ellipsoidalVincenty as V, trf

    t = Tests(__file__, __version__, trf)
    t.testTrf(K.Cartesian, K.LatLon, 'Karney')
    t.testTrf(N.Cartesian, N.LatLon, 'Nvector')
    t.testTrf(V.Cartesian, V.LatLon, 'Vincenty')
    t.testEpoch()
    t.results()
    t.exit()
