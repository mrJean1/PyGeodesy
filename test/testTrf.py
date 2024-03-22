
# -*- coding: utf-8 -*-

# Test L{trf} I{Terrestrial Reference Frame} implementation.

# All tests transcoded from Chris Veness (C) 2006-2022 U{Geodesy tools for conversions between
# reference frames<https://www.Movable-Type.co.UK/scripts/geodesy-library.html>} JavaScript original
# <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-ellipsoidal-referenceframe-tests.js>

__all__ = ('Tests',)
__version__ = '24.03.12'

from bases import GeodSolve, isPython2, TestsBase

from pygeodesy import date2epoch, Epoch, epoch2date, F_D, F_DMS, fstr, RefFrames, \
                      TRFError, trfTransform0, trfTransforms, trfXform, Vector3d


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
        self.test('convertRefFrame', x, '[3980574.395, -102.214, 4966829.941]')

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
        self.test('Geodetic', x.toStr(form=F_DMS, prec=5), '23°40′12.44582″S, 133°53′07.84795″E, +603.34m', nl=1)  # +603.3361m
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

        # Tutorial <https://EPNCB.OMA.Be/_productsservices/coord_trans>
        def _ts(x, c, r1, r2, e, **exhaust):
            for i, T in enumerate(trfTransforms(r1, r2, epoch=e, **exhaust)):
                x = T.toRefFrame(x)
                d = c - x
                d = max(map(abs, d.xyz))
                t = '%s max %.5g, epoched %s' % (T.name, d, fstr(T.Xform.epoched, prec=3))
                self.test('transform%s' % (i,), t, t)

        c = Cartesian(4027894.006, 307045.600, 4919474.910)
        x = c.toRefFrame(RefFrames.ITRF2005, RefFrames.ITRF91, 2007)
        self.test('EUREF EX1', x.toStr(prec=4), '[4027894.0444, 307045.6209, 4919474.8613]', known=isPython2, nl=1)
        _ts(x, c, 'ITRF2005', 'ITRF91', 2007, exhaust=True)

        x = x.toRefFrame(RefFrames.ITRF91, RefFrames.ITRF2005, 2007)
        self.test('EUREF EX2', x.toStr(prec=-4), c.toStr(prec=-4), known=isPython2)
        _ts(x, c, 'ITRF91', 'ITRF2005', 2007)

        x = c.toRefFrame(RefFrames.ETRF2000, RefFrames.ITRF2000, 2012)
        self.test('EUREF EX4', x.toStr(prec=-4), '[4027894.3559, 307045.2508, 4919474.6447]', known=isPython2)
        _ts(x, c, 'ETRF2000', 'ITRF2000', 2012)

        x = c.toRefFrame(RefFrames.ITRF2014, RefFrames.ETRF2000, 2012)
        self.test('EUREF EX5', x.toStr(prec=-4), '[4027894.3662, 307045.2530, 4919474.6263]', known=isPython2 or abs(x.x - 4027894.3662) < 0.75)
        _ts(x, c, 'ITRF2014', 'ETRF2000', 2012)

        # Altamimi, Z. U{"EUREF Technical Note 1: Relationship and Transformation between the International and
        # the European Terrestrial Reference Systems"<https://ERTS89.ENSG,IFN.Fr/pub/EUREF-TN-1.pdf>} Appendix B.
        c = Cartesian(4027893.6719, 307045.9064, 4919475.1704, reframe=RefFrames.ITRF2014, epoch=2010.0)
        x = c.toRefFrame(RefFrames.ETRF2014, epoch=2010)
        self.test('Case 1A', x.toStr(prec=4), '[4027893.9619, 307045.5481, 4919474.9553]', nl=1)  # was .9620, .5480
        c = Cartesian(4027893.6812, 307045.9082, 4919475.1547, reframe=RefFrames.ITRF2000)
        x = c.toRefFrame(RefFrames.ETRF2000, epoch=2010)
        self.test('Case 1B', x.toStr(prec=4), '[4027894.0054, 307045.5938, 4919474.9083]')  # was .0053, .5939
        c = Cartesian(4027893.5358, 307046.0740, 4919475.2748, reframe=RefFrames.ITRF2014, epoch=2020.0)
        x = c.toRefFrame(RefFrames.ETRF2014, epoch=2020)
        self.test('Case 2A', x.toStr(prec=4), '[4027893.9639, 307045.545, 4919474.9573]')
        c = Cartesian(4027893.5505, 307046.0772, 4919475.245, reframe=RefFrames.ITRF2000)
        x = c.toRefFrame(RefFrames.ETRF2000, epoch=2020)
        self.test('Case 2B', x.toStr(prec=4), '[4027894.0033, 307045.5889, 4919474.9041]')  # was .0033, .5889, .9047

        try:
            t = LatLon(0, 0, reframe='ITRF2000')
        except TypeError as x:
            t = str(x)
        self.test('TypeError', t, "type(reframe) ('ITRF2000'): not a RefFrame", nl=1)

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
            h = trfTransform0(c.reframe, RefFrames.ITRF2014, epoch=c.epoch, epoch2=2018.8)
            self.test('TransformXform', h.name, h.name)
            e = x - X
            self.test('Delta (m)', e.toStr(prec=6), '[0.01, 0.01, 0.01]', known=max(map(abs, e.xyz)) < 0.5)
            e = e.length
            self.test('Error (m)', e, '0.01', prec=6, known=e < 1.5)
            r = x.epoch - c.epoch
            self.test('Epoch range', r, '14.0', prec=3, known=True)
            r = x.toRefFrame(c.reframe, epoch2=c.epoch)
            self.test(_n(r), r.toStr(prec=-6), t, known=(r - c).length < 1)
            h = trfTransform0(x.reframe, c.reframe, epoch=x.epoch, epoch2=c.epoch)
            self.test('TransformXform', h.name, h.name)

        def _t(x, dX):
            return '%s@%s %s %s' % (x.reframe.name, x.reframe.epoch, x.epoch, dX.toStr(prec=8))

        x = RefFrames.ITRF2014.toRefFrame(X, RefFrames.ITRF2020, epoch=2018.8, epoch2=2024.31)
        self.test('toRefFrame1', _t(x, x - X), "ITRF2020@2015 2024.31 [0.0031474, 0.00210534, -0.00125667]", nl=1)
        x = RefFrames.ITRF2020.Xform('ITRF2014').toRefFrame(X, epoch=2018.8, epoch2=2024.32)
        self.test('toRefFrame2', _t(x, X - x), "ITRF2014@2010 2024.32 [0.0031474, 0.00210634, -0.00125867]")  # flipped
        T = trfTransform0('ITRF2014', RefFrames.ITRF2020)  # get Transform from reframe, no epoch or epoch2 ...
        self.test('transform0', T.name, T.name)
        t = repr(T.Xform)
        self.test('transform0X', t, t)
        x = T.toRefFrame(X)
        self.test('toRefFrame3', _t(x, x - X), "ITRF2020@2015 2010 [0.0031474, 0.00067434, 0.00160533]")
        x = T.transform(X.x, X.y, X.z)
        self.test('transform2x', x.toStr(prec=-6), '(4160476.488147, 653193.021674, 4774604.781605)')
        v = T.velocities()
        self.test('transform2v', v.toStr(prec=-3), '(0.004, 0.003, 0.004)', known=True)

        def _t4(c, x, p, X=None):
            x = Vector3d(x)
            e = c - x
            d = max(map(abs, e.xyz))
            X = (', epoched %s' % (fstr(X.epoched, prec=3),)) if X else ''
            return c.toStr(prec=p), x.toStr(prec=p), d, '%s max %.4g%s' % (e.toStr(prec=9), d, X)

        # Alamimi, Z. Example 1 <http://ETRS89.ENSG.IGN.FR/pub/EUREF-TN-1-Jan-31-2024.pdf>
        f, v, r1 =        (4027893.6750, 307045.9069, 4919475.1721), (-0.01361,  0.01686,  0.01024), 'ITRF2020'
        for t, w, r2 in (((4027893.9585, 307045.5550, 4919474.9619), ( 0.00011,  0.00011,  0.00024), 'ETRF2020'),
                         ((4027893.6719, 307045.9064, 4919475.1704), (-0.01361,  0.01676,  0.01044), 'ITRF2014'),
                         ((4027893.9620, 307045.5480, 4919474.9553), ( 0.00020, -0.00030,  0.00020), 'ETRF2014'),
                         ((4027893.6812, 307045.9082, 4919475.1547), (-0.01307,  0.01690,  0.00908), 'ITRF2000'),
                         ((4027894.0053, 307045.5939, 4919474.9083), (-0.00020, -0.00050, -0.00036), 'ETRF2000')):
            T = trfTransform0(r1, r2, epoch=2010)
            x = T.toStr(prec=9)
            self.test('transform0', x, x, nl=1)
            c, v = T.transform2(*(f + v), Vector=Vector3d)
            c, x, d, e = _t4(c, t, -4, T.Xform)
            self.test('transform2c', c, x, known=d < 1.e-4)
            self.test('    Error2c', e, e)
            v, x, d, e = _t4(v, w, -5)
            self.test('transform2v', v, x, known=d < 1.0 or d > 100.0)
            self.test('    Error2v', e, e)
            v = T.velocities(Vector=Vector3d)
            v, x, d, e = _t4(v, w, -5)
            self.test('transform0v', v, x, known=d < 1.0)
            self.test('    Error0v', e, e)
            f, r1, v = t, r2, w

        f,    r1 = (4027893.9620, 307045.5480, 4919474.9553,    0.00020, -0.00030,  0.00020), 'ETRF2014'
        t, w, r2 = (4027893.6812, 307045.9082, 4919475.1547), (-0.01307,  0.01690,  0.00908), 'ITRF2000'
        for i, T in enumerate(trfTransforms(r1, r2, epoch=2010)):
            n = 'transform' + str(i)
            x = T.toStr(prec=9)
            self.test(n + '/', x, x, nl=1)
            c, v = T.transform2(*f, Vector=Vector3d)
            c, x, d, e = _t4(c, t, -4, T.Xform)
            self.test('transform2c', c, x, known=d < 0.2)
            self.test('    Error2c', e, e)
            v, x, d, e = _t4(v, w, -5)
            self.test('transform2v', v, x, known=d < 1.0 or d > 100.0)
            self.test('    Error2v', e, e)

        # Alamimi, Z. Example 2 <http://ETRS89.ENSG.IGN.FR/pub/EUREF-TN-1-Jan-31-2024.pdf>
        f, r1 =        (4027893.5389, 307046.0755, 4919475.2745), 'ITRF2020'
        for t, r2 in (((4027893.9574, 307045.5561, 4919474.9643), 'ETRF2020'),
                      ((4027893.5358, 307046.0740, 4919475.2748), 'ITRF2014'),
                      ((4027893.9639, 307045.5450, 4919474.9573), 'ETRF2014'),
                      ((4027893.5505, 307046.0772, 4919475.2456), 'ITRF2000'),
                      ((4027894.0033, 307045.5889, 4919474.9047), 'ETRF2000')):
            T = trfTransform0(r1, r2, epoch=2020)
            x = T.toStr(prec=9)
            self.test('transform0', x, x, nl=1)
            c, v = T.transform2(*f, Vector=Vector3d)
            c, x, d, e = _t4(c, t, -4, T.Xform)
            self.test('transform2c', c, x, known=d < 1.e-4)
            self.test('    Error2c', e, e)
            v = v.toStr(prec=-5)
            self.test('transform2v', v, v)
            f, r1 = t, r2

        # courtesy GGaessler <https://GitHub.com/mrJean1/PyGeodesy/issues/80>
        t, w = (4160476.485, 653193.021, 4774604.780), (0.004, 0.003, 0.004)
        for f, r1, e in (((4160476.944, 653192.600, 4774604.455), 'ETRF89',   1989),
                         ((4160476.952, 653192.582, 4774604.441), 'ETRF2000', 2000),
                         ((4160476.674, 653192.806, 4774604.648), 'ETRF2008', 2005)):
            for i, T in enumerate(trfTransforms(r1, 'ITRF2014', epoch=e, epoch2=2018.8)):
                n = 'transform' + str(i)
                x = T.toStr(prec=9)
                self.test(n + '*', x, x, nl=1)
                c, v = T.transform2(*(f + w), Vector=Vector3d)
                c, x, d, e = _t4(c, t, -4, T.Xform)
                self.test('transform2c', c, x, known=d < 0.2)
                self.test('    Error2c', e, e)
                v, x, d, e = _t4(v, w, -5)
                self.test('transform2v', v, x, known=d < 0.03)
                self.test('    Error2v', e, e)

        t = str(T)
        self.test('toTransform', t, t, nl=1)
        t = str(T.toTransform(epoch1=T.Xform.epoch))
        self.test('toTransform', t, t)

        # courtesy GGaessler <https://GitHub.com/mrJean1/PyGeodesy/issues/80>
        p = LatLon('48 46 36.89676N', '8 55 21.25713E', height=476.049, reframe=RefFrames.ETRF89, epoch=1989)
        self.test('Issue80', p.toStr(form='D', prec=8), '48.77691577°N, 008.92257142°E, +476.05m', nl=1)
        x = p.toRefFrame(RefFrames.ITRF2014, epoch2=2018.8)
        self.test('Issue80', x.toStr(form='D', prec=8), '48.77692147°N, 008.92257868°E, +476.09m')
        t = x.toRefFrame(RefFrames.ETRF89, epoch=1989)  # 48.7769157667°N, 008.9225714250°E, +476.049m
        d = p.vincentysTo(t)
        self.test('Issue80', t.toStr(prec=6), '48°46′36.915134″N, 008°55′21.285094″E, +476.10', known=d < 0.9)
        self.test('Issue80', d, '0.01', prec=3, known=d < 0.9)

        c = p.toCartesian(Cartesian=Cartesian)  # reframe=RefFrames.ETRF89, epoch=1989
        self.test('Issue80', c.toStr(prec=6), '[4160476.944064, 653192.600457, 4774604.455385]')
        x = c.toRefFrame(RefFrames.ITRF2014, epoch2=2018.8)
        self.test('Issue80', x.toStr(prec=6), '[4160476.415504, 653193.057171, 4774604.903316]')
        t = x.toRefFrame(RefFrames.ETRF89, epoch=1989)
        self.test('Issue80', t.toStr(prec=6), '[4160476.467901, 653193.103647, 4774604.868888]')
        p = t.toLatLon(LatLon=LatLon)  # reframe=RefFrames.ETRF89, epoch=1989
        self.test('Issue80', p.toStr(prec=6), '48°46′36.915133″N, 008°55′21.285094″E, +476.10m')

        p = LatLon('48 46 36.91314N', '8 55 21.28095E', height=476.049, reframe=RefFrames.ITRF2014, epoch=2018.8)
        self.test('Issue80', p.toStr(form='D', prec=8), '48.77692032°N, 008.92257804°E, +476.05m', nl=1)
        x = p.toRefFrame(RefFrames.ETRF89, epoch2=1989)
        self.test('Issue80', x.toStr(form='D', prec=8, nl=1), '48.77691971°N, 008.92257856°E, +476.06m')
        t = x.toRefFrame(RefFrames.ITRF2014, epoch=2018.8)   # 48.7769203167°N, 008.9225780417°E, +476.037m
        d = p.vincentysTo(t)
        self.test('Issue80', t.toStr(prec=4), '48°46′36.9131″N, 008°55′21.28095″E, +476.05m', known=d < 0.9)
        self.test('Issue80', d, '0.01', prec=3, known=d < 0.9)

        c = p.toCartesian(Cartesian=Cartesian)  # reframe=RefFrames.ITRF2014, epoch=2018.8
        self.test('Issue80', c.toStr(prec=6), '[4160476.492633, 653193.021888, 4774604.78885]')
        x = c.toRefFrame(RefFrames.ETRF89, epoch2=1989)
        self.test('Issue80', x.toStr(prec=6), '[4160476.54503, 653193.068365, 4774604.754422]')
        t = x.toRefFrame(RefFrames.ITRF2014, epoch=2018.8)
        self.test('Issue80', t.toStr(prec=6), '[4160476.016469, 653193.525079, 4774605.202353]')
        t = t.toLatLon(LatLon=LatLon)  # reframe=RefFrames.ITRF2014, epoch=2018.8
        self.test('Issue80', t.toStr(prec=4), '48°46′36.9315″N, 008°55′21.3089″E, +476.10m')

        # <https://Geodesy.NOAA.gov/TOOLS/Htdp/Using_HTDP.pdf> equations (5)-(11)
        X = trfXform('ITRF96', 'NAD83')
        self.test(X.name, X.toRepr(), repr(RefFrames.ITRF96.Xform('NAD83')), nl=1)
        h = X.toHelmert()
        self.test(X.name, h, "name='ITRF96@1997xNAD83', tx=0.991, ty=-0.19072, tz=-0.5129, s1=1.0,"
                                                      " rx=1.2503e-07, ry=4.6785e-08, rz=5.6529e-08, s=0.0,"
                                                      " sx=0.02579, sy=0.00965, sz=0.01166")
        t = X.toTransform(2007)  # to show the 10-year deltas for rx, ry and rz
        self.test(X.name, t, "name='ITRF96@1997xNAD83@2007', tx=0.991, ty=-0.19072, tz=-0.5129, s1=1.0,"
                                                           " rx=1.2761e-07, ry=1.0797e-08, rz=5.4997e-08, s=0.0,"
                                                           " sx=0.026322, sy=0.002227, sz=0.011344")

        def _d(h, t):
            return '%.5g' % ((h - t) * 0.1)  # per year, see (8)-(10)

        t = 'rx=%s, ry=%s, rz=%s' % (_d(h.rx, t.rx), _d(h.ry, t.ry), _d(h.rz, t.rz))
        self.test(X.name, t, "rx=-2.5792e-10, ry=3.5988e-09, rz=1.532e-10")

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
        self.test('repr', repr(r), t, nt=1)

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
        self.test(e.toRepr() + '.std_repr', e.std_repr, False, nl=1)
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
