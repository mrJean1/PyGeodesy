
# -*- coding: utf-8 -*-

'''NavLab examples

This page illustrates implementations of the examples from
<https://www.NavLab.net/nvector>.  Values are chosen to match
those used in nvector.readthedocs.org.  Tests marked with
# +++ are additional, not present in the original examples.
'''
__all__ = ()
__version__ = '21.05.17'

if __name__ == '__main__':

    from base import GeodSolve, geographiclib, TestsBase

    from pygeodesy import Datums, F_D, ellipsoidalExact, \
                          ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    class Examples(TestsBase):  # overload test()
        def test(self, ex, name, *args, **kwds):
            name = 'Example %s %s' % (ex, name)
            TestsBase.test(self, name, *args, **kwds)

    def destination(m, x):
        a = m.LatLon(80, -90)  # +++
        b = a.destination(1000, 200)
        n = '%s(%s)' % (destination.__name__, m.__name__)
        t.test(8, n, b.toStr(F_D), x)

    t = Examples(__file__, __version__)

# Example 1: A and B to delta
    a = ellipsoidalNvector.LatLon(1, 2, 3)  # defaults to WGS-84
    b = ellipsoidalNvector.LatLon(4, 5, 6)
    delta = a.deltaTo(b)
    t.test(1, 'delta', delta, '[331730.863, 332998.501, 17398.304]')
    t.test(1, 'delta', delta.toRepr(prec=3), '[L:470357.384, B:45.109°, E:-2.12°]')  # DEPRECATED
    t.test(1, 'elevation', delta.elevation, -2.1198, fmt='%.4f')
    t.test(1, 'bearing', delta.bearing, 45.109, fmt='%.3f')  # 45.109°
    t.test(1, 'length', delta.length, 470357.384, fmt='%.3f')  # 470357.384 m

# Example 2: B and delta to C*
    n = ellipsoidalNvector.Nvector(1, 2, 3, 400, Datums.WGS72)
#   t.test(2, 'Nvector', n.toStr(prec=3), '[1.0, 2.0, 3.0, +400.00]')
    b = n.toLatLon()
    t.test(2, 'LatLon', b.toStr(F_D, prec=3), '53.301°N, 063.435°E, +400.00m')
    t.test(2, 'toNvector', b.toNvector().toStr(prec=3), '(0.267, 0.535, 0.802, +400.00)')
    delta = ellipsoidalNvector.Ned(3000, 2000, 100)
    t.test(2, 'delta', delta, '[3000.0, 2000.0, 100.0]')  # ++
    t.test(2, 'delta', delta.toRepr(prec=3), '[L:3606.938, B:33.69°, E:-1.589°]')  # DEPRECATED
    c = b.destinationNed(delta)  # JSname: destinationPoint
    t.test(2, 'destinationNed', c.toStr(F_D), '53.327726°N, 063.464965°E, +299.138m', known=True)

    a = ellipsoidalNvector.LatLon(49.66618, 3.45063)  # ++
    b = ellipsoidalNvector.LatLon(48.88667, 2.37472)  # ++
    delta = a.deltaTo(b)  # ++
    t.test(2, 'delta', delta.toStr(prec=0), '[-86126, -78900, 1069]')  # ++
    t.test(2, 'delta', delta.toRepr(prec=3), '[L:116807.681, B:222.493°, E:-0.524°]')  # DEPRECATED
    c = a.destinationNed(delta)  # JSname: destinationPoint, c.height = -9.3e-10
    t.test(2, 'destinationNed', c.toStr(F_D), '48.88667°N, 002.37472°E', known=True)

# Example 3: ECEF-vector to geodetic latitude
    c = ellipsoidalNvector.Cartesian(0.9*6371e3, -1.0*6371e3, 1.1*6371e3)
#   t.test(3, 'Cartesian', c, '[5733900.0, -6371000.0, 7008100.0]')
    p = c.toLatLon()
    t.test(3, 'toLatLon', p.toStr(F_D, prec=3), '39.379°N, 048.013°W, +4702059.83m')

# Example 4: Geodetic latitude to ECEF-vector
    p = ellipsoidalNvector.LatLon(1, 2, 3)
    c = p.toCartesian()
    t.test(4, 'toCartesian', c.toStr(prec=3), '[6373290.277, 222560.201, 110568.827]')

# Example 5: Surface distance
    a = sphericalNvector.LatLon(88, 0)
    b = sphericalNvector.LatLon(89, -170)
    dist = a.distanceTo(b)
    t.test(5, 'distanceTo', dist, 332457, fmt='%.0f')  # 332,457 m == 332.5 km

# Example 6: Interpolated position
    a = sphericalNvector.LatLon(89, 0)
    b = sphericalNvector.LatLon(89, 180)
    p = a.intermediateChordTo(b, 0.6)
    t.test(6, 'intermediateChordTo', p.toStr(F_D), '89.799981°N, 180.0°E')
    p = a.intermediateTo(b, 0.6)
    t.test(6, 'intermediateTo', p.toStr(F_D), '89.8°N, 180.0°E')

    a = sphericalNvector.LatLon(52.205, 0.119)  # +++
    b = sphericalNvector.LatLon(48.857, 2.351)
    p = a.intermediateChordTo(b, 0.25)
    t.test(6, 'intermediateChordTo', p.toStr(F_D), '51.372294°N, 000.707192°E')  # 51.3723°N, 000.7072°E
    p = a.intermediateTo(b, 0.25)
    t.test(6, 'intermediateTo', p.toStr(F_D), '51.372084°N, 000.707337°E')

# Example 7: Mean position
    points = [sphericalNvector.LatLon(90,   0),
              sphericalNvector.LatLon(60,  10),
              sphericalNvector.LatLon(50, -20)]
    mean = sphericalNvector.meanOf(points)  # XXX meanOf
    t.test(7, 'meanOf', mean.toStr(F_D, prec=4), '67.2362°N, 006.9175°W')
#   t.test(7, 'meanOfLatLon', mean.__class__, "<class 'sphericalNvector.LatLon'>")  # ++

# Example 8: A and azimuth/distance to B
    destination(sphericalNvector,      '79.991549°N, 090.017698°W')
    destination(sphericalTrigonometry, '79.991549°N, 090.017698°W')
    destination(ellipsoidalVincenty,   '79.991584°N, 090.017621°W')
    if geographiclib:
        from pygeodesy import ellipsoidalKarney
        destination(ellipsoidalKarney, '79.991584°N, 090.017621°W')
    destination(ellipsoidalExact,      '79.991584°N, 090.017621°W')
    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        destination(ellipsoidalGeodSolve, '79.991584°N, 090.017621°W')

# Example 9: Intersection of two paths
    a1 = sphericalNvector.LatLon(10, 20)
    a2 = sphericalNvector.LatLon(30, 40)
    b1 = sphericalNvector.LatLon(50, 60)
    b2 = sphericalNvector.LatLon(70, 80)
    c = sphericalNvector.intersection(a1, a2, b1, b2)
    t.test(9, 'intersection', c, '40.318643°N, 055.901868°E')

# Example 10: Cross track distance
    a1 = sphericalNvector.LatLon( 0, 0)
    a2 = sphericalNvector.LatLon(10, 0)
    b = sphericalNvector.LatLon(1, 0.1)
    c = b.crossTrackDistanceTo(a1, a2)
    t.test(10, 'crossTrackDistance', c, 11118, fmt='%.0f')  # 11,118 m == 11.12 km

# <https://GitHub.com/chrisveness/geodesy/blob/master/latlon-nvector-ellipsoidal.js>
    d = ellipsoidalNvector.toNed(116809.178, 222.493, -0.5416)
    TestsBase.test(t, 'toNed', d.toStr(prec=1), '[-78901.1, -86126.6, 1104.1]')  # [N:-86126.6, E:-78901.1, D:1104.1]'
    TestsBase.test(t, 'bearing',   d.bearing, '227.507',  fmt='%.3f')  # '222.493'
    TestsBase.test(t, 'elevation', d.elevation, '-0.5416', fmt='%.4f')
    TestsBase.test(t, 'length',    d.length, '116809.178',  fmt='%.3f')
    v = d.toVector3d()
    TestsBase.test(t, 'toVector3d', v.toStr(prec=1), '(-86126.6, -78901.1, -1104.1)')  # 1104.1

    t.results()
    t.exit()
