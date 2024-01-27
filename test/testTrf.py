
# -*- coding: utf-8 -*-

# Test L{trf} I{Terrestrial Reference Frame} implementation.

# All tests transcoded from Chris Veness (C) 2006-2022 U{Geodesy tools for conversions between
# reference frames<https://www.Movable-Type.co.UK/scripts/geodesy-library.html>} JavaScript original
# <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-ellipsoidal-referenceframe-tests.js>

__all__ = ('Tests',)
__version__ = '24.01.25'

from bases import GeodSolve, TestsBase

from pygeodesy import date2epoch, Epoch, epoch2date, F_D, F_DMS, RefFrames, TRFError


class Tests(TestsBase):

    def testTrf(self, ellipsoidal_):  # MCCABE 17

        self.subtitle(ellipsoidal_, 'Trf')

        Cartesian = ellipsoidal_.Cartesian
        LatLon    = ellipsoidal_.LatLon

        p = LatLon(51.47788, -0.00147, reframe=RefFrames.ITRF2000)
        x = p.convertRefFrame(RefFrames.ETRF2000)
        self.test('convertRefFrame', x.toStr(form=F_D, prec=8), '51.47787826°N, 000.00147125°W', known=True)

        c = Cartesian(4027893.924, 307041.993, 4919474.294)
        x = c.toLatLon()
        x.reframe = RefFrames.ITRF2000
        self.testCopy(x.reframe)
        self.test('toLatLon', x.toStr(form=F_D, prec=4), '50.7978°N, 004.3592°E, +148.96m')
        c = Cartesian(3980574.247, -102.127, 4966830.065)
        x = c.toRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2000)
        self.test('convertRefFrame', x, '[3980574.395, -102.209, 4966829.941]')  # [3980574.395, -102.214, 4966829.941]

        p = LatLon(0, 0, reframe=RefFrames.ITRF2000)
        x = p.toRefFrame(RefFrames.ITRF2000)
        self.test('Nil', x is p, True)
        self.testCopy(x.reframe)
        c = Cartesian(1, 2, 3)
        x = c.toRefFrame(RefFrames.ITRF2000, RefFrames.ITRF2000)
        self.test('Nil', x is c, True)

        p = LatLon(0, 0, reframe=RefFrames.NAD83)
        x = p.toRefFrame(RefFrames.ITRF2014)  # # via ITRF2000
        self.test('reframe', x.reframe == RefFrames.ITRF2014, True)
        x = p.toRefFrame(RefFrames.NAD83)
        self.test('Roundtrip', x == p, True)
        self.test('reframe', x.reframe == RefFrames.NAD83, True)
        self.testCopy(x.reframe)

        # Dawson, J. & Woods, A., Appendix A, Journal of Applied Geodesy 4 (2010)
        p = LatLon('23°40′12.41482″S', '133°53′7.86712″E', height=603.2562, reframe=RefFrames.ITRF2005, epoch=2010.4559)
        x = p.toRefFrame(RefFrames.GDA94)
        self.test('Geodetic', x.toStr(form=F_DMS, prec=5), '23°40′12.44582″S, 133°53′07.84795″E, +603.34m')  # +603.3361m
        c = x.toCartesian()
        self.test('Cartesian', c.toStr(prec=4), '[-4052051.7614, 4212836.1945, -2545106.0147]')

        x = p.toRefFrame(RefFrames.GDA94)  # epoch 2010.4559
        self.test('Geodetic', x.toStr(form=F_DMS, prec=5), '23°40′12.44582″S, 133°53′07.84795″E, +603.34m')  # +603.3361m'
        c = x.toCartesian()
        self.test('Cartesian', c.toStr(prec=4), '[-4052051.7614, 4212836.1945, -2545106.0147]')

        x = x.toRefFrame(RefFrames.ITRF2005)  # epoch 2010.4559
        self.test('Roundtrip', x.toStr(form=F_DMS, prec=5), '23°40′12.41482″S, 133°53′07.86712″E, +603.26m')  # +603.2562m

        # <https://GitHub.com/OSGeo/proj.4/blob/2aaf53/test/gie/more_builtins.gie#L357>
        c = Cartesian(3370658.37800, 711877.31400, 5349787.08600)  # Proj4 Onsala observatory
        x = c.toRefFrame(RefFrames.ITRF93, RefFrames.ITRF2000, 2017)
        self.test('GNSStrans', x.toStr(prec=5), '[3370658.18892, 711877.42369, 5349787.1243]')  # accurate to within 0.02 mm

        # <https://www.NGS.NOAA.gov/cgi-bin/ds_mark.prl?PidBox=kg0640>
        p = LatLon('39 13 26.71220', '098 32 31.74540', height=573.961, reframe=RefFrames.NAD83, epoch=2010.0)
        c = p.toCartesian()  # NGS Data Sheet Meades Ranch
        self.test('Cartesian', c, '[-734972.563, 4893188.492, 4011982.811]')

        # <https://EPNCB.OMA.Be/_productsservices/coord_trans> (tutorial)
        c = Cartesian(4027894.006, 307045.600, 4919474.910)
        x = c.toRefFrame(RefFrames.ITRF91, RefFrames.ITRF2005, 2007)
        self.test('EUREF C1', x.toStr(prec=4), '[4027894.0444, 307045.6209, 4919474.8613]')
        x = c.toRefFrame(RefFrames.ITRF91, RefFrames.ITRF2005, 2007)
        self.test('EUREF C2', x.toStr(prec=4), '[4027894.0444, 307045.6209, 4919474.8613]')
        x = c.toRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2000, 2012)
        self.test('EUREF C4', x.toStr(prec=4), '[4027894.356, 307045.2501, 4919474.6447]')  # [4027894.3559, 307045.2508, 4919474.6447]
        x = c.toRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2014, 2012)
        self.test('EUREF C5', x.toStr(prec=4), '[4027894.3662, 307045.253, 4919474.6263]')

        # Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
        # the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Appendix B.
        c = Cartesian(4027893.6719, 307045.9064, 4919475.1704, reframe=RefFrames.ITRF2014, epoch=2010.0)
        x = c.toRefFrame(RefFrames.ETRF2014, epoch=2010)
        self.test('Case 1A', x.toStr(prec=4), '[4027893.9619, 307045.5481, 4919474.9553]')  # was .9620, .5480
        c = Cartesian(4027893.6812, 307045.9082, 4919475.1547, reframe=RefFrames.ITRF2000)
        x = c.toRefFrame(RefFrames.ETRF2000, epoch=2010)
        self.test('Case 1B', x.toStr(prec=4), '[4027894.0054, 307045.5938, 4919474.9083]')  # was .0053, .5939
        c = Cartesian(4027893.5358, 307046.0740, 4919475.2748, reframe=RefFrames.ITRF2014, epoch=2020.0)
        x = c.toRefFrame(RefFrames.ETRF2014, epoch=2020)
        self.test('Case 2A', x.toStr(prec=4), '[4027893.9639, 307045.545, 4919474.9573]')  # was .5450
        c = Cartesian(4027893.5505, 307046.0772, 4919475.245, reframe=RefFrames.ITRF2000)
        x = c.toRefFrame(RefFrames.ETRF2000, epoch=2020)
        self.test('Case 2B', x.toStr(prec=4), '[4027894.0036, 307045.585, 4919474.9041]')  # was .0033, .5889, .9047

        try:
            t = LatLon(0, 0, reframe='ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe) ('ITRF2000'): not a RefFrame")

        try:
            t = LatLon(0, 0, reframe=RefFrames.ITRF2000, epoch=1899)
        except TRFError as x:
            t = str(x)
        self.test(TRFError.__name__, t, "epoch (1899): below 1900.0 limit")

        try:
            t = LatLon(0, 0, reframe=RefFrames.ITRF2000).convertRefFrame('ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe2) ('ITRF2000'): not a RefFrame")

        try:
            t = LatLon(0, 0).toRefFrame(RefFrames.ITRF2000)
        except TRFError as x:
            t = str(x)
        self.test('TRFError', t, 'no conversion: LatLon(00°00′00.0″N, 000°00′00.0″E).reframe MISSING')

        c = Cartesian(0, 0, 0)
        try:
            t = c.toRefFrame('ITRF2000', RefFrames.ITRF2000)
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe2) ('ITRF2000'): not a RefFrame")

        try:
            t = c.toRefFrame(RefFrames.ITRF2000, 'ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe) ('ITRF2000'): not a RefFrame")

#       try:
#           t = c.convertRefFrame(RefFrames.ITRF2000, RefFrames.ITRF2000, '2000')
#       except TypeError as x:
#           t = str(x)
#       self.test('TypeError', t, "type(epoch) ('2000'): not scalar")

        def _n(c):
            return '%s@%s' % (c.reframe.name, c.epoch)

        X = Cartesian(4160476.485, 653193.021, 4774604.780)  # expected
        for c in (Cartesian(4160476.944, 653192.600, 4774604.455, reframe=RefFrames.ETRF89,   epoch=1989),
                  Cartesian(4160476.952, 653192.582, 4774604.441, reframe=RefFrames.ETRF2000, epoch=2000),
                  Cartesian(4160476.674, 653192.806, 4774604.648, reframe=RefFrames.ITRF2008, epoch=2005)):
            t = c.toStr(prec=-6)
            self.test(_n(c), t, t, known=True, nl=1)
            x = c.toRefFrame(RefFrames.ITRF2014, epoch2=2018.8)
            self.test(_n(x), x.toStr(prec=-3), X, known=True)
            e = x - X
            self.test('Delta (m)', e.toStr(prec=6), '[0.01, 0.01, 0.01]', known=max(map(abs, e.xyz)) < 0.5)
            e = e.length
            self.test('Error (m)', e, '0.01', prec=6, known=e < 1.5)
            r = x.epoch - c.epoch
            self.test('Epoch range', r, '14.0', prec=3, known=True)
            x = x.toRefFrame(c.reframe, epoch2=c.epoch)
            self.test(_n(x), x.toStr(prec=-6), t)

        x = RefFrames.ITRF2014.toRefFrame(X, RefFrames.ITRF2020, epoch=2018.8, epoch2=2024.1)  # coverage
        x = (x - X).toStr(prec=6)
        self.test('toRefFrame', x, x, nl=1)

        # courtesy GGaessler <https://github.com/mrJean1/PyGeodesy/issues/80>
        p = LatLon('48 46 36.89676N', '8 55 21.25713E', height=476.049, reframe=RefFrames.ETRF89, epoch=1989)
        self.test('Issue80', p.toStr(form='D', prec=10), '48.7769157667°N, 008.922571425°E, +476.05m', nl=1)
        x = p.toRefFrame(RefFrames.ITRF2014, epoch2=2018.8)
        self.test('Issue80', x.toStr(form='D', prec=10), '48.7769169572°N, 008.9225709519°E, +476.09m')
        t = x.toRefFrame(RefFrames.ETRF89, epoch=1989)  # 48.7769157667°N, 008.9225714250°E, +476.049m
        self.test('Issue80', t.toStr(prec=6), '48°46′36.898876″N, 008°55′21.257278″E, +476.10m')

        c = p.toCartesian(Cartesian=Cartesian)  # reframe=RefFrames.ETRF89, epoch=1989
        self.test('Issue80', c.toStr(prec=6), '[4160476.944064, 653192.600457, 4774604.455385]')
        x = c.toRefFrame(RefFrames.ITRF2014, epoch2=2018.8)
        self.test('Issue80', x.toStr(prec=6), '[4160476.875697, 653192.554523, 4774604.571079]')
        t = x.toRefFrame(RefFrames.ETRF89, epoch=1989)
        self.test('Issue80', t.toStr(prec=6), '[4160476.928094, 653192.601, 4774604.536651]')
        p = t.toLatLon(LatLon=LatLon)  # reframe=RefFrames.ETRF89, epoch=1989
        self.test('Issue80', p.toStr(prec=6), '48°46′36.898876″N, 008°55′21.257278″E, +476.10m')

        p = LatLon('48 46 36.91314N', '8 55 21.28095E', height=476.049, reframe=RefFrames.ITRF2014, epoch=2018.8)
        self.test('Issue80', p.toStr(form='D', prec=10), '48.7769203167°N, 008.9225780417°E, +476.05m', nl=1)
        x = p.toRefFrame(RefFrames.ETRF89, epoch2=1989)
        self.test('Issue80', x.toStr(form='D', prec=10, nl=1), '48.7769191262°N, 008.9225785148°E, +476.01m')
        t = x.toRefFrame(RefFrames.ITRF2014, epoch=2018.8)    # 48.7769203167°N, 008.9225780417°E, +476.037m
        self.test('Issue80', t.toStr(prec=6), '48°46′36.929398″N, 008°55′21.308766″E, +476.05m')

        c = p.toCartesian(Cartesian=Cartesian)  # reframe=RefFrames.ITRF2014, epoch=2018.8
        self.test('Issue80', c.toStr(prec=6), '[4160476.492633, 653193.021888, 4774604.78885]')
        x = c.toRefFrame(RefFrames.ETRF89, epoch2=1989)
        self.test('Issue80', x.toStr(prec=6), '[4160476.561, 653193.067822, 4774604.673156]')
        t = x.toRefFrame(RefFrames.ITRF2014, epoch=2018.8)
        self.test('Issue80', t.toStr(prec=6), '[4160476.03244, 653193.524536, 4774605.121087]')
        p = t.toLatLon(LatLon=LatLon)  # reframe=RefFrames.ITRF2014, epoch=2018.8
        self.test('Issue80', p.toStr(prec=6), '48°46′36.929398″N, 008°55′21.308766″E, +476.05m')

    def testEpoch(self):

        try:  # coverage
            e = date2epoch(None, 1, 2)
            self.test('epoch', e, TRFError.__name__, nl=1)
        except TRFError as x:
            t = str(x)
            self.test('TRFError', t, t, nl=1)

        r = RefFrames.GDA94
        t = r.toStr()
        self.test('toStr', t, "name='GDA94', epoch=1994, datum=Datums.GRS80")
        self.test('str', str(r),t)
        t = r.toRepr()
        self.test('toStr2', t, "RefFrame(name='GDA94', epoch=1994, datum=Datums.GRS80)")
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

    from pygeodesy import ellipsoidalExact, ellipsoidalKarney, \
                          ellipsoidalNvector, ellipsoidalVincenty, trf

    t = Tests(__file__, __version__, trf)
    t.testTrf(ellipsoidalNvector)
    t.testTrf(ellipsoidalVincenty)
    t.testTrf(ellipsoidalKarney)
    t.testTrf(ellipsoidalExact)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testTrf(ellipsoidalGeodSolve)

    t.testEpoch()
    t.results()
    t.exit()
