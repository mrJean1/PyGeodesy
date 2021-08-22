
# -*- coding: utf-8 -*-

u'''Ellipsoidal and spherical earth models.

Classes L{a_f2Tuple}, L{Ellipsoid} and L{Ellipsoid2}, an L{Ellipsoids} registry and
a dozen functions to convert I{equatorial} radius, I{polar} radius, I{eccentricities},
I{flattenings} and I{inverse flattening}.

See module L{datums} for more information and other details.

@var Ellipsoids.Airy1830: Ellipsoid(name='Airy1830', a=6377563.396, b=6356256.90923729, f_=299.3249646, f=0.00334085, f2=0.00335205, n=0.00167322, e=0.08167337, e2=0.00667054, e22=0.00671533, e32=0.00334643, A=6366914.60892522, L=10001126.0807165, R1=6370461.23374576, R2=6370459.65470808, R3=6370453.30994572, Rbiaxial=6366919.065224, Rtriaxial=6372243.45317691)
@var Ellipsoids.AiryModified: Ellipsoid(name='AiryModified', a=6377340.189, b=6356034.44793853, f_=299.3249646, f=0.00334085, f2=0.00335205, n=0.00167322, e=0.08167337, e2=0.00667054, e22=0.00671533, e32=0.00334643, A=6366691.77461988, L=10000776.05340819, R1=6370238.27531284, R2=6370236.69633043, R3=6370230.35179013, Rbiaxial=6366696.2307627, Rtriaxial=6372020.43236847)
@var Ellipsoids.Australia1966: Ellipsoid(name='Australia1966', a=6378160, b=6356774.71919531, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, A=6367471.84853228, L=10002001.39064442, R1=6371031.5730651, R2=6371029.9824858, R3=6371023.59124344, Rbiaxial=6367476.337459, Rtriaxial=6372820.40754721)
@var Ellipsoids.Bessel1841: Ellipsoid(name='Bessel1841', a=6377397.155, b=6356078.962818, f_=299.1528128, f=0.00334277, f2=0.00335398, n=0.00167418, e=0.08169683, e2=0.00667437, e22=0.00671922, e32=0.00334836, A=6366742.52023395, L=10000855.76443237, R1=6370291.09093933, R2=6370289.51012659, R3=6370283.15821523, Rbiaxial=6366746.98155108, Rtriaxial=6372074.29334012)
@var Ellipsoids.CPM1799: Ellipsoid(name='CPM1799', a=6375738.7, b=6356671.92557493, f_=334.39, f=0.00299052, f2=0.00299949, n=0.0014975, e=0.07727934, e2=0.0059721, e22=0.00600798, e32=0.00299499, A=6366208.88184784, L=10000017.52721564, R1=6369383.10852498, R2=6369381.8434158, R3=6369376.76247022, Rbiaxial=6366212.45090321, Rtriaxial=6370977.3559758)
@var Ellipsoids.Clarke1866: Ellipsoid(name='Clarke1866', a=6378206.4, b=6356583.8, f_=294.97869821, f=0.00339008, f2=0.00340161, n=0.00169792, e=0.08227185, e2=0.00676866, e22=0.00681478, e32=0.00339582, A=6367399.68916978, L=10001888.04298286, R1=6370998.86666667, R2=6370997.240633, R3=6370990.70659881, Rbiaxial=6367404.2783313, Rtriaxial=6372807.62791066)
@var Ellipsoids.Clarke1880: Ellipsoid(name='Clarke1880', a=6378249.145, b=6356514.86954978, f_=293.465, f=0.00340756, f2=0.00341921, n=0.00170669, e=0.0824834, e2=0.00680351, e22=0.00685012, e32=0.00341337, A=6367386.64398051, L=10001867.55164747, R1=6371004.38651659, R2=6371002.74366963, R3=6370996.1419165, Rbiaxial=6367391.2806777, Rtriaxial=6372822.52526083)
@var Ellipsoids.Clarke1880IGN: Ellipsoid(name='Clarke1880IGN', a=6378249.2, b=6356515, f_=293.46602129, f=0.00340755, f2=0.0034192, n=0.00170668, e=0.08248326, e2=0.00680349, e22=0.00685009, e32=0.00341336, A=6367386.73667336, L=10001867.69724907, R1=6371004.46666667, R2=6371002.82383112, R3=6370996.22212395, Rbiaxial=6367391.37333829, Rtriaxial=6372822.59907505)
@var Ellipsoids.Clarke1880Mod: Ellipsoid(name='Clarke1880Mod', a=6378249.145, b=6356514.96582849, f_=293.4663, f=0.00340755, f2=0.0034192, n=0.00170668, e=0.08248322, e2=0.00680348, e22=0.00685009, e32=0.00341335, A=6367386.69207875, L=10001867.62720001, R1=6371004.4186095, R2=6371002.77577708, R3=6370996.17408252, Rbiaxial=6367391.32873482, Rtriaxial=6372822.54926891)
@var Ellipsoids.Delambre1810: Ellipsoid(name='Delambre1810', a=6376428, b=6355957.92616372, f_=311.5, f=0.00321027, f2=0.00322061, n=0.00160772, e=0.08006397, e2=0.00641024, e22=0.0064516, e32=0.00321543, A=6366197.07684334, L=9999998.98395793, R1=6369604.64205457, R2=6369603.18419749, R3=6369597.32739068, Rbiaxial=6366201.19059818, Rtriaxial=6371316.64722284)
@var Ellipsoids.Engelis1985: Ellipsoid(name='Engelis1985', a=6378136.05, b=6356751.32272154, f_=298.2566, f=0.00335282, f2=0.0033641, n=0.00167922, e=0.08181928, e2=0.00669439, e22=0.00673951, e32=0.00335844, A=6367448.17507971, L=10001964.20447208, R1=6371007.80757385, R2=6371006.21707085, R3=6370999.82613573, Rbiaxial=6367452.66379074, Rtriaxial=6372796.59560563)
@var Ellipsoids.Everest1969: Ellipsoid(name='Everest1969', a=6377295.664, b=6356094.667915, f_=300.8017, f=0.00332445, f2=0.00333554, n=0.00166499, e=0.08147298, e2=0.00663785, e22=0.0066822, e32=0.00332998, A=6366699.57839501, L=10000788.3115495, R1=6370228.665305, R2=6370227.10178537, R3=6370220.81951618, Rbiaxial=6366703.99082487, Rtriaxial=6372002.02812501)
@var Ellipsoids.Fisher1968: Ellipsoid(name='Fisher1968', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367463.65604381, L=10001988.52191361, R1=6371022.77908146, R2=6371021.18903735, R3=6371014.79995035, Rbiaxial=6367468.14345752, Rtriaxial=6372811.30979281)
@var Ellipsoids.GEM10C: Ellipsoid(name='GEM10C', a=6378137, b=6356752.31424783, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367449.14582474, L=10001965.7293148, R1=6371008.77141594, R2=6371007.18091936, R3=6371000.79001005, Rbiaxial=6367453.63451765, Rtriaxial=6372797.55596006)
@var Ellipsoids.GRS67: Ellipsoid(name='GRS67', a=6378160, b=6356774.51609071, f_=298.24716743, f=0.00335292, f2=0.0033642, n=0.00167928, e=0.08182057, e2=0.00669461, e22=0.00673973, e32=0.00335854, A=6367471.74706533, L=10002001.2312605, R1=6371031.50536357, R2=6371029.91475409, R3=6371023.52339015, Rbiaxial=6367476.23607738, Rtriaxial=6372820.3568989)
@var Ellipsoids.GRS80: Ellipsoid(name='GRS80', a=6378137, b=6356752.31414035, f_=298.2572221, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367449.14577104, L=10001965.72923046, R1=6371008.77138012, R2=6371007.18088351, R3=6371000.78997414, Rbiaxial=6367453.634464, Rtriaxial=6372797.55593326)
@var Ellipsoids.Helmert1906: Ellipsoid(name='Helmert1906', a=6378200, b=6356818.16962789, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367513.57227074, L=10002066.93013953, R1=6371072.7232093, R2=6371071.13315272, R3=6371064.74401563, Rbiaxial=6367518.05971963, Rtriaxial=6372861.26794141)
@var Ellipsoids.IERS1989: Ellipsoid(name='IERS1989', a=6378136, b=6356751.30156878, f_=298.257, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181922, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367448.13949125, L=10001964.14856985, R1=6371007.76718959, R2=6371006.17669088, R3=6370999.78577297, Rbiaxial=6367452.62819019, Rtriaxial=6372796.55279934)
@var Ellipsoids.IERS1992TOPEX: Ellipsoid(name='IERS1992TOPEX', a=6378136.3, b=6356751.61659215, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367448.44699641, L=10001964.63159783, R1=6371008.07219738, R2=6371006.48170097, R3=6371000.09079236, Rbiaxial=6367452.93568883, Rtriaxial=6372796.85654541)
@var Ellipsoids.IERS2003: Ellipsoid(name='IERS2003', a=6378136.6, b=6356751.85797165, f_=298.25642, f=0.00335282, f2=0.0033641, n=0.00167922, e=0.0818193, e2=0.0066944, e22=0.00673951, e32=0.00335844, A=6367448.71771058, L=10001965.05683465, R1=6371008.35265722, R2=6371006.76215217, R3=6371000.37120877, Rbiaxial=6367453.20642742, Rtriaxial=6372797.14192686)
@var Ellipsoids.Intl1924: Ellipsoid(name='Intl1924', a=6378388, b=6356911.94612795, f_=297, f=0.003367, f2=0.00337838, n=0.00168634, e=0.08199189, e2=0.00672267, e22=0.00676817, e32=0.00337267, A=6367654.50005758, L=10002288.29898944, R1=6371229.31537598, R2=6371227.71133444, R3=6371221.26587487, Rbiaxial=6367659.02704315, Rtriaxial=6373025.77129687)
@var Ellipsoids.Intl1967: Ellipsoid(name='Intl1967', a=6378157.5, b=6356772.2, f_=298.24961539, f=0.0033529, f2=0.00336418, n=0.00167926, e=0.08182023, e2=0.00669455, e22=0.00673967, e32=0.00335852, A=6367469.33894446, L=10001997.44859308, R1=6371029.06666667, R2=6371027.47608389, R3=6371021.08482752, Rbiaxial=6367473.827881, Rtriaxial=6372817.9027631)
@var Ellipsoids.Krassovski1940: Ellipsoid(name='Krassovski1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367558.49687498, L=10002137.49754285, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.69367439, Rbiaxial=6367562.98435553, Rtriaxial=6372906.23027515)
@var Ellipsoids.Krassowsky1940: Ellipsoid(name='Krassowsky1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367558.49687498, L=10002137.49754285, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.69367439, Rbiaxial=6367562.98435553, Rtriaxial=6372906.23027515)
@var Ellipsoids.Maupertuis1738: Ellipsoid(name='Maupertuis1738', a=6397300, b=6363806.28272251, f_=191, f=0.0052356, f2=0.00526316, n=0.00262467, e=0.10219488, e2=0.01044379, e22=0.01055402, e32=0.00524931, A=6380564.13011837, L=10022566.69846922, R1=6386135.42757417, R2=6386131.54144847, R3=6386115.8862823, Rbiaxial=6380575.11882818, Rtriaxial=6388943.03218495)
@var Ellipsoids.Mercury1960: Ellipsoid(name='Mercury1960', a=6378166, b=6356784.28360711, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367479.62923643, L=10002013.61254591, R1=6371038.76120237, R2=6371037.17115427, R3=6371030.78205124, Rbiaxial=6367484.1166614, Rtriaxial=6372827.29640037)
@var Ellipsoids.Mercury1968Mod: Ellipsoid(name='Mercury1968Mod', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367463.65604381, L=10001988.52191361, R1=6371022.77908146, R2=6371021.18903735, R3=6371014.79995035, Rbiaxial=6367468.14345752, Rtriaxial=6372811.30979281)
@var Ellipsoids.NWL1965: Ellipsoid(name='NWL1965', a=6378145, b=6356759.76948868, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, A=6367456.87366841, L=10001977.86818326, R1=6371016.58982956, R2=6371014.999254, R3=6371008.60802667, Rbiaxial=6367461.36258457, Rtriaxial=6372805.42010473)
@var Ellipsoids.OSU86F: Ellipsoid(name='OSU86F', a=6378136.2, b=6356751.51693008, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367448.3471653, L=10001964.47478349, R1=6371007.97231003, R2=6371006.38181364, R3=6370999.99090513, Rbiaxial=6367452.83585765, Rtriaxial=6372796.75662978)
@var Ellipsoids.OSU91A: Ellipsoid(name='OSU91A', a=6378136.3, b=6356751.6165948, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367448.44699773, L=10001964.6315999, R1=6371008.07219827, R2=6371006.48170186, R3=6371000.09079324, Rbiaxial=6367452.93569015, Rtriaxial=6372796.85654607)
@var Ellipsoids.Plessis1817: Ellipsoid(name='Plessis1817', a=6376523, b=6355862.93325557, f_=308.64, f=0.00324002, f2=0.00325055, n=0.00162264, e=0.08043347, e2=0.00646954, e22=0.00651167, e32=0.00324527, A=6366197.15710739, L=9999999.11003639, R1=6369636.31108519, R2=6369634.82608583, R3=6369628.85999668, Rbiaxial=6366201.34758009, Rtriaxial=6371364.26393357)
@var Ellipsoids.SGS85: Ellipsoid(name='SGS85', a=6378136, b=6356751.30156878, f_=298.257, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181922, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367448.13949125, L=10001964.14856985, R1=6371007.76718959, R2=6371006.17669087, R3=6370999.78577297, Rbiaxial=6367452.62819019, Rtriaxial=6372796.55279934)
@var Ellipsoids.SoAmerican1969: Ellipsoid(name='SoAmerican1969', a=6378160, b=6356774.71919531, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, A=6367471.84853228, L=10002001.39064442, R1=6371031.5730651, R2=6371029.98248581, R3=6371023.59124344, Rbiaxial=6367476.337459, Rtriaxial=6372820.40754721)
@var Ellipsoids.Sphere: Ellipsoid(name='Sphere', a=6371008.771415, b=6371008.771415, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, A=6371008.771415, L=10007557.17611675, R1=6371008.771415, R2=6371008.771415, R3=6371008.771415, Rbiaxial=6371008.771415, Rtriaxial=6371008.771415)
@var Ellipsoids.SphereAuthalic: Ellipsoid(name='SphereAuthalic', a=6371000, b=6371000, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, A=6371000, L=10007543.39801029, R1=6371000, R2=6371000, R3=6371000, Rbiaxial=6371000, Rtriaxial=6371000)
@var Ellipsoids.SpherePopular: Ellipsoid(name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, A=6378137, L=10018754.17139462, R1=6378137, R2=6378137, R3=6378137, Rbiaxial=6378137, Rtriaxial=6378137)
@var Ellipsoids.Struve1860: Ellipsoid(name='Struve1860', a=6378298.3, b=6356657.14266956, f_=294.73, f=0.00339294, f2=0.00340449, n=0.00169935, e=0.0823065, e2=0.00677436, e22=0.00682056, e32=0.00339869, A=6367482.31832549, L=10002017.83655714, R1=6371084.58088985, R2=6371082.95208988, R3=6371076.40691418, Rbiaxial=6367486.91530791, Rtriaxial=6372894.90029454)
@var Ellipsoids.WGS60: Ellipsoid(name='WGS60', a=6378165, b=6356783.28695944, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, A=6367478.63091189, L=10002012.0443814, R1=6371037.76231981, R2=6371036.17227197, R3=6371029.78316994, Rbiaxial=6367483.11833616, Rtriaxial=6372826.29723739)
@var Ellipsoids.WGS66: Ellipsoid(name='WGS66', a=6378145, b=6356759.76948868, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, A=6367456.87366841, L=10001977.86818326, R1=6371016.58982956, R2=6371014.999254, R3=6371008.60802667, Rbiaxial=6367461.36258457, Rtriaxial=6372805.42010473)
@var Ellipsoids.WGS72: Ellipsoid(name='WGS72', a=6378135, b=6356750.52001609, f_=298.26, f=0.00335278, f2=0.00336406, n=0.0016792, e=0.08181881, e2=0.00669432, e22=0.00673943, e32=0.0033584, A=6367447.24862383, L=10001962.74919858, R1=6371006.84000536, R2=6371005.24953886, R3=6370998.8587507, Rbiaxial=6367451.7372317, Rtriaxial=6372795.60727472)
@var Ellipsoids.WGS84: Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, A=6367449.14582341, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000916, Rbiaxial=6367453.63451633, Rtriaxial=6372797.5559594)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import copysign0, isfinite, isint, _xinstanceof
from pygeodesy.errors import _AssertionError, _ValueError
from pygeodesy.fmath import cbrt, cbrt2, fdot, fhorner, fpowers, Fsum, fsum_, \
                            hypot, hypot_, hypot1, hypot2, sqrt3
from pygeodesy.interns import EPS, EPS0, EPS02, EPS1, INF, NN, PI4, PI_2, R_M, _a_, \
                             _Airy1830_, _AiryModified_, _Bessel1841_, _Clarke1866_, \
                             _Clarke1880IGN_, _DOT_, _1_EPS, _EPStol as _TOL, _f_, \
                             _finite_, _float as _F, _floatuple as _T, _GRS80_, _height_, \
                             _Intl1924_, _Krassovski1940_, _Krassowsky1940_, _lat_, \
                             _meridional_, _negative_, _not_, _null_, _prime_vertical_, \
                             _radius_, _Sphere_, _SPACE_, _vs_, _WGS72_, _WGS84_, \
                             _0_0, _0_5, _1_0, _2_0, _4_0, _90_0
from pygeodesy.interns import _0_25, _3_0  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _lazyNamedEnumItem as _lazy, _NamedEnum, \
                            _NamedEnumItem, _NamedTuple, _Pass
from pygeodesy.namedTuples import Distance2Tuple, Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_Property_RO, Property_RO, property_doc_
from pygeodesy.streprs import Fmt, fstr, instr, strs, unstr
from pygeodesy.units import Bearing_, Distance, Float, Float_, Height, Lam_, Lat, Meter, \
                            Meter2, Meter3, Phi, Phi_, Radius, Radius_, Scalar
from pygeodesy.utily import atand, atan2b, atan2d, degrees90, m2km, m2NM, m2SM, \
                            m2radians, radians2m, sincos2d

from math import asinh, atan, atanh, cos, degrees, exp, radians, sin, sinh, sqrt, tan

R_M  = Radius(R_M =R_M)            # mean (spherical) earth radius (C{meter})
R_MA = Radius(R_MA=_F(6378137.0))  # equatorial earth radius (C{meter}), WGS84, EPSG:3785
R_MB = Radius(R_MB=_F(6356752.3))  # polar earth radius (C{meter}), WGS84, EPSG:3785
R_KM = Radius(R_KM=_F(m2km(R_M)))  # mean (spherical) earth radius (C{KM}, kilo meter)
R_NM = Radius(R_NM=_F(m2NM(R_M)))  # mean (spherical) earth radius (C{NM}, nautical miles)
R_SM = Radius(R_SM=_F(m2SM(R_M)))  # mean (spherical) earth radius (C{SM}, statute miles)
# See <https://www.EdWilliams.org/avform.htm>,
# <https://www.DTIC.mil/dtic/tr/fulltext/u2/a216843.pdf> and
# <https://GitHub.com/NASA/MultiDop/blob/master/src/share/man/man3/geog_lib.3>
# based on International Standard Nautical Mile of 1,852 meter (1' latitude)
R_FM = Radius(R_FM=_F(6371000.0))        # former FAI Sphere earth radius (C{meter})
R_GM = Radius(R_GM=_F(6371230.0))        # Avg. radius, distance to geoid surface (C{meter})
R_VM = Radius(R_VM=_F(6366707.0194937))  # Aviation/Navigation earth radius (C{meter})
# R_ = Radius(R_  =_F(6372797.560856))   # XXX some other earth radius???

__all__ = _ALL_LAZY.ellipsoids
__version__ = '21.08.21'

_f_0_0   = Float(f =_0_0)
_f__0_0  = Float(f_=_0_0)
# like U{WGS84_f()<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Constants.html>}
_f_WGS84 = Float(f = 1 / (1000000000 / 298257223563))


def _aux(lat, inverse, auxLat, clip=90):
    '''Return named auxiliary latitude in C{degrees}.
    '''
    return Lat(lat, clip=clip, name=_lat_ if inverse else auxLat)


def _4Ecef(this, Ecef):
    '''Return an ECEF converter.
    '''
    from pygeodesy.ecef import EcefKarney, EcefVeness, EcefYou

    if Ecef is None:
        Ecef = EcefKarney
    else:
        _xinstanceof(EcefKarney, EcefVeness, EcefYou, Ecef=Ecef)
    return Ecef(this, name=this.name)  # datum or ellipsoid


def _s2_c2(phi):
    '''(INTERNAL) Return 2-tuple C{(sin(B{phi})**2, cos(B{phi})**2)}.
    '''
    if phi:
        s2 = sin(phi)**2
        if s2 > 0:
            c2 = _1_0 - s2
            if c2 > 0:
                if c2 < _1_0:
                    return s2, c2
            else:
                return _1_0, _0_0
    return _0_0, _1_0


class a_f2Tuple(_NamedTuple):
    '''2-Tuple C{(a, f)} specifying an ellipsoid by I{equatorial}
       radius C{a} in C{meter} and scalar I{flattening} C{f}.

       @see: Class L{Ellipsoid2}.
    '''
    _Names_ = (_a_,   _f_)  # name 'f' not 'f_'
    _Units_ = (_Pass, _Pass)

    def __new__(cls, a, f, **name):
        '''New L{a_f2Tuple} ellipsoid specification.

           @arg a: Equatorial radius (C{scalar} > 0).
           @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

           @return: An L{a_f2Tuple}C{(a, f)} instance.

           @raise UnitError: Invalid B{C{a}} or B{C{f}}.

           @note: C{abs(B{f}) < EPS} is forced to C{B{f}=0}, I{spherical}.
                  Negative C{B{f}} produces a I{prolate} ellipsoid.
        '''
        a = Radius_(a=a)
        f = Float_( f=f, low=None, high=EPS1)
        if abs(f) < EPS:  # force spherical
            f = _f_0_0
        return _NamedTuple.__new__(cls, a, f, **name)

    @Property_RO
    def b(self):
        '''Get the I{polar} radius (C{meter}), M{a * (1 - f)}.
        '''
        return a_f2b(self.a, self.f)  # PYCHOK .a and .f

    @Property_RO
    def f_(self):
        '''Get the I{inverse} flattening (C{float}), M{1 / f} == M{a / (a - b)}.
        '''
        return f2f_(self.f)  # PYCHOK .f


class Circle4Tuple(_NamedTuple):
    '''4-Tuple C{(radius, height, lat, beta)} of the C{radius} and C{height},
       both conventionally in C{meter} of a parallel I{circle of latitude} at
       (geodetic) latitude C{lat} and the I{parametric (or reduced) auxiliary
       latitude} C{beta}, both in C{degrees90}.

       The C{height} is the (signed) distance along the z-axis between the
       parallel and the equator.  At near-polar C{lat}s, the C{radius} is C{0},
       the C{height} is the ellipsoid's (signed) polar radius and C{beta}
       equals C{lat}.
    '''
    _Names_ = (_radius_, _height_, _lat_, 'beta')
    _Units_ = ( Radius,   Height,   Lat,   Lat)


class Curvature2Tuple(_NamedTuple):
    '''2-Tuple C{(meridional, prime_vertical)} of radii of curvature,
       both in C{meter}, conventionally.
    '''
    _Names_ = (_meridional_, _prime_vertical_)
    _Units_ = ( Meter,        Meter)


class Ellipsoid(_NamedEnumItem):
    '''Ellipsoid with I{equatorial} and I{polar} radii, I{flattening}, I{inverse
       flattening} and other, often used, I{cached} attributes, supporting
       I{oblate} and I{prolate} ellipsoidal and I{spherical} earth models.
    '''
    _a  = 0  # equatorial radius, semi-axis (C{meter})
    _b  = 0  # polar radius, semi-axis (C{meter}): a * (f - 1) / f
    _f  = 0  # (1st) flattening: (a - b) / a
    _f_ = 0  # inverse flattening: 1 / f = a / (a - b)

    _KsOrder = 8     # Krüger series order (4, 6 or 8)
    _Math    = None  # cached karney._wrapped.Math module

    def __init__(self, a, b=None, f_=None, name=NN):
        '''New L{Ellipsoid} from I{equatorial} and I{polar} radius or
           I{equatorial} radius and I{inverse flattening}.

           @arg a: Equatorial radius, semi-axis (C{meter}).
           @arg b: Optional, polar radius, semi-axis (C{meter}).
           @arg f_: Inverse flattening: M{a / (a - b)} (C{float} >>> 1.0).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}}, B{C{b}} or B{C{f_}}.

           @note: M{abs(f_) > 1 / EPS} or M{abs(1 / f_) < EPS} is forced
                  to M{1 / f_ = 0}, spherical.
        '''
        try:
            a = Radius_(a=a)  # low=EPS
            if not isfinite(a):
                raise ValueError(_not_(_finite_))

            if b:
                b  = Radius_(b=b)  # low=EPS
                f  = a_b2f(a, b)
                f_ = f2f_(f) if f_ is None else Float(f_=f_)
            elif f_:
                f_ = Float(f_=f_)
                b  = a_f_2b(a, f_)  # a * (f_ - 1) / f_
                f  = a_b2f(a, b)
            else:  # only a, spherical
                f = f_ = 0
                b = a  # superfluous

            if not isfinite(b):
                raise ValueError(_not_(_finite_))

            if abs(f) < EPS or a == b or not f_:  # spherical
                b  =  a
                f  = _f_0_0
                f_ = _f__0_0
            elif f > EPS1:  # sanity check, see .ecef.Ecef.__init__
                raise _ValueError(f=f)

        except (TypeError, ValueError) as x:
            t = instr(self, a=a, b=b, f_=f_, name=name)
            raise _ValueError(t, txt=str(x))

        self._a  = a
        self._b  = b
        self._f  = f
        self._f_ = f_

        self._register(Ellipsoids, name)

        if f and f_:  # see .test/testEllipsoidal.py
            self._assert(_1_0 / f,  f_=f_, eps=_TOL)
            self._assert(_1_0 / f_, f =f,  eps=_TOL)
        self._assert(self.b2_a2, e12=self.e12, eps=EPS)

    def __eq__(self, other):
        '''Compare this and an other ellipsoid.

           @arg other: The other ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Ellipsoid) and
                                  self.a == other.a and
                                 (self.b == other.b or self.f == other.f))

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self._a

    equatoradius = a  # = Requatorial

    @Property_RO
    def a2(self):
        '''Get the I{equatorial} radius I{squared} (C{meter**2}), M{a**2}.
        '''
        return Meter2(a2=self.a**2)

    @Property_RO
    def a2_(self):
        '''Get the inverse of the I{equatorial} radius I{squared} (C{meter**2}), M{1 / a**2}.
        '''
        return Float(a2_=_1_0 / self.a2)

    @Property_RO
    def a_b(self):
        '''Get the ratio I{equatorial} over I{polar} radius (C{float}), M{a / b} == M{1 / (1 - f)}.
        '''
        return Float(a_b=self.a / self.b if self.f else _1_0)

    @Property_RO
    def a2_b(self):
        '''Get the I{polar} meridional radius of curvature (C{meter}), M{a**2 / b}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
                 and U{Moritz, H. (1980), Geodetic Reference System 1980
                 <https://WikiPedia.org/wiki/Earth_radius#cite_note-Moritz-2>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of curvature
                  <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}, see L{c2} and
                  L{R2} or L{Rauthalic}.
        '''
        return Radius(a2_b=self.a2 / self.b if self.f else self.a)  # = rocPolar

    @Property_RO
    def a2_b2(self):
        '''Get the ratio I{equatorial} over I{polar} radius I{squared} (C{float}), M{(a / b)**2}
           == M{1 / (1 - e**2)} == M{1 / (1 - e2)} == M{1 / e12}.
        '''
        return Float(a2_b2=self.a_b**2 if self.f else _1_0)

    @Property_RO
    def a_f(self):
        '''Get the I{equatorial} radius and I{flattening} (L{a_f2Tuple}).
        '''
        return a_f2Tuple(self.a, self.f, name=self.name)

    @Property_RO
    def A(self):
        '''Get the UTM I{meridional (or rectifying)} radius (C{meter}).
        '''
        A, n = self.a, self.n
        if n:
            n1 = _1_0 + n
            if n1:  # use 6 n**2 terms, half-way between the _KsOrder's 4, 6, 8
                # <https://GeographicLib.SourceForge.io/html/tmseries30.html>
                # <https://GeographicLib.SourceForge.io/html/transversemercator.html> and
                # <https://www.MyGeodesy.id.AU/documents/Karney-Krueger%20equations.pdf> (3)
                A = Radius(A=A / n1 * (fhorner(n**2, 1048576, 262144, 16384, 4096, 1600, 784, 441) / 1048576))
        return A

    @Property_RO
    def _albersCyl(self):
        '''(INTERNAL) Helper for C{auxAuthalic}.
        '''
        from pygeodesy.albers import AlbersEqualAreaCylindrical
        return AlbersEqualAreaCylindrical(datum=self, name=self.name)

    @Property_RO
    def AlphaKs(self):
        '''Get the I{Krüger} U{Alpha series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}C{-tuple}).
        '''
        return self._Kseries(  # XXX int/int quotients may require  from __future__ import division
            #   n    n**2   n**3      n**4         n**5            n**6                 n**7                     n**8
            _T(1/2, -2/3,   5/16,    41/180,    -127/288,       7891/37800,         72161/387072,        -18975107/50803200),
                 _T(13/48, -3/5,    557/1440,    281/630,   -1983433/1935360,       13769/28800,         148003883/174182400),     # PYCHOK unaligned
                        _T(61/240, -103/140,   15061/26880,   167603/181440,    -67102379/29030400,       79682431/79833600),      # PYCHOK unaligned
                               _T(49561/161280, -179/168,    6601661/7257600,       97445/49896,      -40176129013/7664025600),    # PYCHOK unaligned
                                            _T(34729/80640, -3418889/1995840,    14644087/9123840,      2605413599/622702080),     # PYCHOK unaligned
                                                        _T(212378941/319334400, -30705481/10378368,   175214326799/58118860800),   # PYCHOK unaligned
                                                                            _T(1522256789/1383782400, -16759934899/3113510400),    # PYCHOK unaligned
                                                                                                  _T(1424729850961/743921418240))  # PYCHOK unaligned

    @Property_RO
    def area(self):
        '''Get the ellipsoid's surface area (C{meter**2}), M{4 * PI * c2}.

           @see: L{areax}, L{c2} and L{R2}.
        '''
        return Meter2(area=self.c2 * PI4)

    @Property_RO
    def areax(self):
        '''Get the ellipsoid's surface area (C{meter**2}), M{4 * PI * c2x},
           more accurate for very I{oblate} ellipsoids.

           @see: L{area}, L{c2x}, L{R2x} and L{GeodesicExact}.
        '''
        return Meter2(areax=self.c2x * PI4)

    def _assert(self, val, eps=_TOL, f0=_0_0, **name_value):
        '''(INTERNAL) Assert a C{name=value} vs C{val}.
        '''
        for n, v in name_value.items():
            if abs(v - val) > eps:
                t = (v, _vs_, val)
                t = _SPACE_.join(strs(t, prec=12, fmt=Fmt.g))
                t =  Fmt.EQUAL(self._DOT_(n), t)
                raise _AssertionError(t, txt=Fmt.exceeds_eps(eps))
            return Float(v if self.f else f0, name=n)
        raise _AssertionError(unstr(self._DOT_(self._assert.__name__), val,
                                    eps=eps, f0=f0, **name_value))

    def auxAuthalic(self, lat, inverse=False):
        '''Compute the I{authalic} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{authalic}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{authalic} and
                           return the geodetic latitude (C{bool}).

           @return: The I{authalic} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/AuthalicLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Authalic latitude
                 <https://WikiPedia.org/wiki/Latitude#Authalic_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, p 16.

        '''
        if self.f:
            f = self._albersCyl._tanf if inverse else self._albersCyl._txif  # PYCHOK attr
            lat = atand(f(tan(Phi_(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxAuthalic.__name__)

    def auxConformal(self, lat, inverse=False):
        '''Compute the I{conformal} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{conformal}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{conformal} and
                           return the geodetic latitude (C{bool}).

           @return: The I{conformal} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ConformalLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Conformal latitude
                 <https://WikiPedia.org/wiki/Latitude#Conformal_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 15-16.
        '''
        if self.f:
            f = self.es_tauf if inverse else self.es_taupf  # PYCHOK attr
            lat = atand(f(tan(Phi_(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxConformal.__name__)

    def auxGeocentric(self, lat, inverse=False):
        '''Compute the I{geocentric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{geocentric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the geocentric and
                           return the I{geocentric} latitude (C{bool}).

           @return: The I{geocentric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/GeocentricLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Geocentric latitude
                 <https://WikiPedia.org/wiki/Latitude#Geocentric_latitude>}, and
                 U{Snyder<<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 17-18.
        '''
        if self.f:
            f = self.a2_b2 if inverse else self.b2_a2
            lat = atand(f * tan(Phi_(lat)))
        return _aux(lat, inverse, Ellipsoid.auxGeocentric.__name__)

    def auxIsometric(self, lat, inverse=False):
        '''Compute the I{isometric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{isometric}) latitude (C{degrees}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{isometric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{isometric} (or geodetic) latitude in C{degrees}.

           @note: The I{isometric} latitude for geodetic C{+/-90} is far
                  outside the C{[-90..+90]} range but the inverse
                  thereof is the original geodetic latitude.

           @see: U{Inverse-/IsometricLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Isometric latitude
                 <https://WikiPedia.org/wiki/Latitude#Isometric_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 15-16.
        '''
        if self.f:
            r = Phi_(lat, clip=0)
            lat = degrees( atan(self.es_tauf(sinh(r))) if inverse else
                          asinh(self.es_taupf(tan(r))))
        # clip=0, since auxIsometric(+/-90) is far outside [-90..+90]
        return _aux(lat, inverse, Ellipsoid.auxIsometric.__name__, clip=0)

    def auxParametric(self, lat, inverse=False):
        '''Compute the I{parametric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{parametric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{parametric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{parametric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ParametricLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Parametric latitude
                 <https://WikiPedia.org/wiki/Latitude#Parametric_(or_reduced)_latitude>},
                 and U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, p 18.
        '''
        if self.f:
            lat = self._beta(Lat(lat), inverse=inverse)
        return _aux(lat, inverse, Ellipsoid.auxParametric.__name__)

    auxReduced = auxParametric  # synonyms

    def auxRectifying(self, lat, inverse=False):
        '''Compute the I{rectifying} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{rectifying}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{rectifying} and
                           return the geodetic latitude (C{bool}).

           @return: The I{rectifying} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/RectifyingLatitude<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Rectifying latitude
                 <https://WikiPedia.org/wiki/Latitude#Rectifying_latitude>}, and
                 U{Snyder<https://Pubs.USGS.gov/pp/1395/report.pdf>}, pp 16-17.
        '''
        if self.f:
            lat = Lat(lat)
            if 0 < abs(lat) < _90_0:
                if inverse:
                    e = self._elliptic_e22
                    lat = degrees90(e.fEinv(e.cE * lat / _90_0))
                    lat = self.auxParametric(lat, inverse=True)
                else:
                    lat = _90_0 * self.Llat(lat) / self.L
        return _aux(lat, inverse, Ellipsoid.auxRectifying.__name__)

    @Property_RO
    def b(self):
        '''Get the I{polar} radius, semi-axis (C{meter}).
        '''
        return self._b

    polaradius = b  # = Rpolar

    @Property_RO
    def b_a(self):
        '''Get the ratio I{polar} over I{equatorial} radius (C{float}), M{b / a} == M{1 - f}.
        '''
        return self._assert(self.b / self.a, b_a=_1_0 - self.f, f0=_1_0)

    @Property_RO
    def b2(self):
        '''Get the I{polar} radius I{squared} (C{float}), M{b**2}.
        '''
        return Meter2(b2=self.b**2)

    @Property_RO
    def b2_a(self):
        '''Get the I{equatorial} meridional radius of curvature (C{meter}), M{b**2 / a}, see C{rocMeridional}C{(0)}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return Radius(b2_a=self.b2 / self.a if self.f else self.b)

    @Property_RO
    def b2_a2(self):
        '''Get the ratio I{polar} over I{equatorial} radius I{squared} (C{float}), M{(b / a)**2}
           == M{(1 - f)**2} == M{1 - e**2} == C{e12}.
        '''
        return Float(b2_a2=self.b_a**2 if self.f else _1_0)

    def _beta(self, lat, inverse=False):
        '''(INTERNAL) Get the I{parametric (or reduced) auxiliary latitude} or inverse thereof.
        '''
        s, c = sincos2d(lat)  # like Karney's tand(lat)
        s *= self.a_b if inverse else self.b_a
        return atan2d(s, c)  # == atand(s / c) if c else copysign0(_90_0, lat)

    @Property_RO
    def BetaKs(self):
        '''Get the I{Krüger} U{Beta series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}C{-tuple}).
        '''
        return self._Kseries(  # XXX int/int quotients may require  from __future__ import division
            #   n    n**2  n**3     n**4        n**5            n**6                 n**7                   n**8
            _T(1/2, -2/3, 37/96,   -1/360,    -81/512,      96199/604800,     -5406467/38707200,      7944359/67737600),
                  _T(1/48, 1/15, -437/1440,    46/105,   -1118711/3870720,       51841/1209600,      24749483/348364800),      # PYCHOK unaligned
                       _T(17/480, -37/840,   -209/4480,      5569/90720,       9261899/58060800,     -6457463/17740800),       # PYCHOK unaligned
                              _T(4397/161280, -11/504,    -830251/7257600,      466511/2494800,     324154477/7664025600),     # PYCHOK unaligned
                                          _T(4583/161280, -108847/3991680,    -8005831/63866880,     22894433/124540416),      # PYCHOK unaligned
                                                      _T(20648693/638668800, -16363163/518918400, -2204645983/12915302400),    # PYCHOK unaligne
                                                                          _T(219941297/5535129600, -497323811/12454041600),    # PYCHOK unaligned
                                                                                              _T(191773887257/3719607091200))  # PYCHOK unaligned

    @deprecated_Property_RO
    def c(self):
        '''DEPRECATED, use property C{R2} or C{Rauthalic}.'''
        return self.R2

    @Property_RO
    def c2(self):
        '''Get the I{authalic} earth radius I{squared} (C{meter**2}).


           @see: L{c2x}, L{area}, L{R2}, L{Rauthalic}, I{Karney's} U{equation 60
                 <https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>} and C++ U{Ellipsoid.Area()
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Ellipsoid.html>},
                 U{Authalic radius<https://WikiPedia.org/wiki/Earth_radius#Authalic_radius>}, U{Surface area
                 <https://WikiPedia.org/wiki/Ellipsoid>} and U{surface area
                 <https://www.Numericana.com/answer/geometry.htm#oblate>}.
        '''
        return self._c2f(False)

    @Property_RO
    def c2x(self):
        '''Get the I{authalic} earth radius I{squared} (C{meter**2}), more accurate for very I{oblate}
           ellipsoids.

           @see: L{c2}, L{areax}, L{R2x}, L{Rauthalicx}, L{GeodesicExact} and I{Karney}'s comments at C++
                 attribute U{GeodesicExact._c2<https://GeographicLib.SourceForge.io/html/GeodesicExact_8cpp_source.html>}.
        '''
        return self._c2f(True)

    def _c2f(self, c2x):
        '''(INTERNAL) Helper for C{.c2} and C{.c2x}.
        '''
        f = self.f
        if f:
            c2, e = self.b2, self.e
            if e > EPS0:
                if f > 0:  # .isOblate
                    c2 *= (asinh(sqrt(self.e22abs)) if c2x else atanh(e)) / e
                elif f < 0:  # .isProlate
                    c2 *= atan(e) / e  # XXX asin?
            c2 = Meter2(c2=(self.a2 + c2) * _0_5)
        else:
            c2 = self.a2
        return c2

    def circle4(self, lat):
        '''Get the equatorial or a parallel I{circle of latitude}.

           @arg lat: Geodetic latitude (C{degrees90}, C{str}).

           @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}
                    instance.

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise TypeError: Invalid B{C{lat}}.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Definition of U{I{p} and I{z} under B{Parametric (or
                 reduced) latitude}<https://WikiPedia.org/wiki/Latitude>}
                 and I{Karney's} C++ U{CircleRadius and CircleHeight
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Ellipsoid.html>}.
        '''
        lat = Lat(lat)
        if lat:
            b = lat
            if abs(lat) < _90_0:
                if self.f:
                    b = self._beta(lat)
                z, r = sincos2d(b)
                r *= self.a
                z *= self.b
            else:  # near-polar
                r, z = _0_0, copysign0(self.b, lat)
        else:  # equator
            r = self.a
            z = lat = b = _0_0
        return Circle4Tuple(r, z, lat, b)

    def degrees2m(self, deg, lat=0):
        '''Convert an angle to the distance along the equator or
           along a parallel of (geodetic) latitude.

           @arg deg: The angle (C{degrees}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Distance (C{meter}, same units as the equatorial
                    and polar radii) or C{0} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{deg}} or B{C{lat}}.
        '''
        return self.radians2m(radians(deg), lat=lat)

    def distance2(self, lat0, lon0, lat1, lon1):
        '''I{Approximate} the distance and (initial) bearing between
           two points based on the U{local, flat earth approximation
           <https://www.EdWilliams.org/avform.htm#flat>} aka U{Hubeny
           <https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           I{Suitable only for distances of several hundred Km or Miles
           and only between points not near-polar}.

           @arg lat0: From latitude (C{degrees}).
           @arg lon0: From longitude (C{degrees}).
           @arg lat1: To latitude (C{degrees}).
           @arg lon1: To longitude (C{degrees}).

           @return: A L{Distance2Tuple}C{(distance, initial)} with C{distance}
                    in same units as this ellipsoid's axes.

           @note: The meridional and prime_vertical radii of curvature are
                  taken and scaled I{at the initial latitude}, see C{roc2}.

           @see: Function L{flatLocal}/L{hubeny}.
        '''
        phi0 = Phi_(lat0=lat0)
        m, n = self.roc2_(phi0, scaled=True)
        m *= Phi_(lat1=lat1) - phi0
        n *= Lam_(lon1=lon1) - Lam_(lon0=lon0)
        return Distance2Tuple(hypot(m, n), atan2b(n, m))

    @Property_RO
    def e(self):
        '''Get the I{(1st) eccentricity} (C{float}), M{sqrt(1 - (b / a)**2))}, see C{a_b2e}.
        '''
        return Float(e=sqrt(self.e2abs) if self.e2 else _0_0)

    eccentricity = e  # eccentricity

    @Property_RO
    def e2(self):
        '''Get the I{(1st) eccentricity squared} (C{float}), M{f * (2 - f)
           == 1 - (b / a)**2}, see C{a_b2e2}.
        '''
        return self._assert(a_b2e2(self.a, self.b), e2=f2e2(self.f))

    eccentricity2 = e2  # eccentricity squared

    @Property_RO
    def e2abs(self):
        '''Get C{abs} value of the I{(1st) eccentricity squared} (C{float}).
        '''
        return abs(self.e2)

    @Property_RO
    def e12(self):
        '''Get 1 less I{(1st) eccentricity squared} (C{float}), M{1 - e**2
           == (1 - f)**2} == M{b**2 / a**2}, see C{b2_a2}.
        '''
        return self._assert((_1_0 - self.f)**2, e12=_1_0 - self.e2, f0=_1_0)  # 1 - e2

    @Property_RO
    def e22(self):
        '''Get the I{2nd eccentricity squared} (C{float}), M{e2 / (1 - e2)
           == e2 / (1 - f)**2 == (a / b)**2 - 1}, see C{a_b2e22}.
        '''
        return self._assert(a_b2e22(self.a, self.b), e22=f2e22(self.f))

    eccentricity2nd2 = e22  # second eccentricity squared

    @Property_RO
    def e22abs(self):
        '''Get C{abs} value of the I{(2nd) eccentricity squared} (C{float}).
        '''
        return abs(self.e22)

    @Property_RO
    def e32(self):
        '''Get the I{3rd eccentricity squared} (C{float}), M{e2 / (2 - e2)
           == (a**2 - b**2) / (a**2 + b**2)}, see C{a_b2e32}.
        '''
        return self._assert(a_b2e32(self.a, self.b), e32=f2e32(self.f))

    eccentricity3rd2 = e32  # third eccentricity squared

    @Property_RO
    def e32abs(self):
        '''Get C{abs} value of the I{(3rd) eccentricity squared} (C{float}).
        '''
        return abs(self.e32)

    @Property_RO
    def e4(self):
        '''Get the I{(1st) eccentricity} to 4th power (C{float}), M{e**4 == e2**2}.
        '''
        return Float(e4=self.e2**2 if self.e2 else _0_0)

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{ellipsoid} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    @Property_RO
    def _elliptic_e22(self):
        '''(INTERNAL) Elliptic helper for C{auxRectifying}, C{L}, C{Llat}.
        '''
        from pygeodesy.elliptic import Elliptic
        return Elliptic(-abs(self.e22))

    def e2s(self, s):
        '''Compute norm M{sqrt(1 - e2 * s**2)}.

           @arg s: Sine value (C{scalar}).

           @return: Norm (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        return sqrt(self.e2s2(s)) if self.e2 else _1_0

    def e2s2(self, s):
        '''Compute M{1 - e2 * s**2}.

           @arg s: S value (C{scalar}).

           @return: Result (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        r = _1_0
        if self.e2:
            try:
                r -= self.e2 * Scalar(s=s)**2
                if r < 0:
                    raise ValueError(_negative_)
            except (TypeError, ValueError) as x:
                t = self._DOT_(Ellipsoid.e2s2.__name__)
                raise _ValueError(t, s, txt=str(x))
        return r

    @Property_RO
    def es(self):
        '''Get the I{(1st) eccentricity signed} (C{float}).
        '''
        # note, self.e is always non-negative
        return Float(es=copysign0(self.e, self.f))  # see .ups

    def es_atanh(self, x):
        '''Compute M{es * atanh(es * x)} where I{es} is the I{signed}
           (1st) eccentricity.

           @raise ValueError: Invalid B{C{x}}.

           @see: Function U{Math::eatanhe<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        # note, self.e is always non-negative
        if self.f > 0:  # .isOblate
            r = atanh(self.e * Scalar(x=x)) * self.e
        elif self.f < 0:  # .isProlate
            r = -atan(self.e * Scalar(x=x)) * self.e
        else:  # .isSpherical
            r = _0_0
        return r

    @Property_RO
    def es_c(self):
        '''Get M{(1 - f) * exp(es_atanh(1))} (C{float}), M{b_a * exp(es_atanh(1))}.
        '''
        return Float(es_c=(exp(self.es_atanh(_1_0)) * self.b_a) if self.f else _1_0)

    def es_tauf(self, taup):
        '''Compute I{Karney}'s U{equations (19), (20) and (21)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: U{Math::tauf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        tp = Scalar(taup=taup)
        if not self.f:  # .isSpherical
            return tp
        tol = max(abs(tp), _1_0) * _TOL
        # To lowest order in e^2, taup = (1 - e^2) * tau, so use starting
        # guess tau = taup / (1 - e^2) = taup * a^2 / b^2.  This starting
        # guess is the geocentric latitude which, to first order in the
        # flattening, is equal to the conformal auxiliary latitude.
        e = self.a2_b2  # == _1_0 / self.e12
        t = tp * e
        T = Fsum(t)
        for _ in range(9):
            a, h = self._es_taupf2(t)
            d = (tp - a) * (e + t**2) / (h * hypot1(a))
            t, d = T.fsum2_(d)
            if abs(d) < tol:
                break
        return t

    def es_taupf(self, tau):
        '''Compute I{Karney}'s U{equations (7), (8) and (9)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: U{Math::taupf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        t = Scalar(tau=tau)
        if self.f:
            t, _ = self._es_taupf2(t)
        return t

    def _es_taupf2(self, t):
        '''(INTERNAL) Return C{(es_taupf(t), hypot1(t))}.
        '''
        h = hypot1(t)
        s = sinh(self.es_atanh(t / h))
        a = hypot1(s) * t - h * s
        return a, h

    @Property_RO
    def f(self):
        '''Get the I{flattening} (C{float}), M{(a - b) / a}, C{0} for spherical, negative for prolate.
        '''
        return self._f

    flattening = f

    @Property_RO
    def f_(self):
        '''Get the I{inverse flattening} (C{float}), M{1 / f} == M{a / (a - b)}, C{0} for spherical, see C{a_b2f_}.
        '''
        return self._f_

    @Property_RO
    def f2(self):
        '''Get the I{2nd flattening} (C{float}), M{(a - b) / b == f / (1 - f)}, C{0} for spherical, see C{a_b2f2}.
        '''
        return self._assert(self.a_b - _1_0, f2=f2f2(self.f))

    flattening2nd = f2

    @Property_RO
    def geodesic(self):
        '''Get this ellipsoid's I{wrapped Karney} U{Geodesic
           <https://GeographicLib.SourceForge.io/html/python/code.html>},
           provided the U{geographiclib<https://PyPI.org/project/geographiclib>}
           package is installed.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        from pygeodesy.karney import _wrapped
        g = _wrapped.Geodesic(self)
        Ellipsoid._Math = _wrapped.Math  # hold
        return g

    @Property_RO
    def _geodesic_Math2(self):
        '''(INTERNAL) Get this ellipsoid's I{wrapped Karney} C{Geodesic}
           and I{Karney}'s C{Math} class, see L{geodesic}.
        '''
        g = self.geodesic
        return g, Ellipsoid._Math

    @Property_RO
    def geodesicx(self):
        '''Get this ellipsoid's L{GeodesicExact}.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        from pygeodesy.geodesicx import GeodesicExact
        return GeodesicExact(self, name=self.name)

    @Property_RO
    def geodsolve(self):
        '''Get this ellipsoid's L{GeodesicSolve}, the I{wrapper} around utility
           U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>},
           provided the path to the C{GeodSolve} executable is specified with env
           variable C{PYGEODESY_GEODSOLVE}.
        '''
        # if not self.isEllipsoidal:
        #     raise _IsnotError(_ellipsoidal_, ellipsoid=self)
        from pygeodesy.geodsolve import GeodesicSolve
        return GeodesicSolve(self, name=self.name)

    def height4(self, xyz, normal=True):
        '''Compute the height of a cartesian above or below and the projection
           on this ellipsoid's surface.

           @arg xyz: The cartesian (C{Cartesian}, L{Ecef9Tuple}, L{Vector3d},
                     L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg normal: If C{True} the projection is the nearest point on the
                          ellipsoid's surface, otherwise the intersection of the
                          radial line to the center and the ellipsoid's surface.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the intersection C{x}, C{y}
                    and C{z} coordinates and height C{h} in C{meter}, conventionally.

           @raise ValueError: Null B{C{xyz}}.

           @raise TypeError: Non-cartesian B{C{xyz}}.

           @see: U{Distance to<https://StackOverflow.com/questions/22959698/distance-from-given-point-to-given-ellipse>}
                 and U{intersection with<https://MathWorld.wolfram.com/Ellipse-LineIntersection.html>} an ellipse.
        '''
        from pygeodesy.vector3d import _otherV3d
        v = _otherV3d(xyz=xyz)
        r =  v.length
        if r < EPS0:  # EPS
            raise _ValueError(xyz=xyz, txt=_null_)

        a, b = self.a, self.b
        if self.isSpherical:
            v = v.times(a / r)
            h = r - a

        elif normal:  # perpendicular to ellipsoid
            x, y = hypot(v.x, v.y), abs(v.z)
            if x < EPS0:  # polar
                z = copysign0(b, v.z)
                v = Vector3Tuple(v.x, v.y, z)
                h = y - b
            elif y < EPS0:  # equatorial
                t = a / r
                v = v.times_(t, t, 0)  # force z=0
                h = x - a
            else:  # normal in 1st quadrant
                x, y = _normal2(x, y, self)
                t, v = v, v.times_(x, x, y)
                h = t.minus(v).length

        else:  # radial to ellipsoid center
            t = hypot_(a * v.z, b * v.x, b * v.y)
            if t < EPS0:  # EPS
                raise _ValueError(xyz=xyz, txt=_null_)
            t = a * b / t
            v = v.times(t)
            h = r * (_1_0 - t)

        return Vector4Tuple(v.x, v.y, v.z, h, name=self.height4.__name__)

    def _hubeny2_(self, phi2, phi1, lam21):
        '''(INTERNAL) like function L{flatLocal_}/L{hubeny_} but
           returning the I{angular} distance in C{radians squared}.
        '''
        m, n = self.roc2_((phi2 + phi1) * _0_5, scaled=True)
        return hypot2(m * (phi2 - phi1), n * lam21) * self.a2_

    @Property_RO
    def isEllipsoidal(self):
        '''Is this model I{ellipsoidal} (C{bool})?
        '''
        return self.f != 0

    @Property_RO
    def isOblate(self):
        '''Is this ellipsoid I{oblate} (C{bool})?  I{Prolate} or
           spherical otherwise.
        '''
        return self.f > 0

    @Property_RO
    def isProlate(self):
        '''Is this ellipsoid I{prolate} (C{bool})?  I{Oblate} or
           spherical otherwise.
        '''
        return self.f < 0

    @Property_RO
    def isSpherical(self):
        '''Is this model I{spherical} (C{bool})?
        '''
        return self.f == 0

    def _Kseries(self, *AB8Ks):
        '''(INTERNAL) Compute the 4-, 6- or 8-th order I{Krüger} Alpha
           or Beta series coefficients per I{Karney}'s U{equations 35
           and 36<https://Arxiv.org/pdf/1002.1417v3.pdf>}.

           @arg AB8Ks: 8-Tuple of 8-th order I{Krüger} Alpha or Beta series
                       coefficient tuples.

           @return: I{Krüger} series coefficients (L{KsOrder}C{-tuple}).

           @see: I{Karney}s 30-th order U{TMseries30
                 <https://GeographicLib.SourceForge.io/html/tmseries30.html>}.
        '''
        k = self.KsOrder
        if self.n:
            ns = fpowers(self.n, k)
            ks = tuple(fdot(AB8Ks[i][:k-i], *ns[i:]) for i in range(k))
        else:
            ks = (_0_0,) * k
        return ks

    @property_doc_(''' the I{Krüger} series' order (C{int}), see properties C{AlphaKs}, C{BetaKs}.''')
    def KsOrder(self):
        '''Get the I{Krüger} series' order (C{int} 4, 6 or 8).
        '''
        return self._KsOrder

    @KsOrder.setter  # PYCHOK setter!
    def KsOrder(self, order):
        '''Set the I{Krüger} series' order.

           @arg order: New I{Krüger} series' order (C{int} 4, 6 or 8).

           @raise ValueError: Invalid B{C{order}}.
        '''
        if not (isint(order) and order in (4, 6, 8)):
            raise _ValueError(order=order)
        if order != self._KsOrder:
            Ellipsoid.AlphaKs._update(self)
            Ellipsoid.BetaKs._update(self)
            self._KsOrder = order

    @Property_RO
    def L(self):
        '''Get the I{quarter meridian} C{L}, aka C{polar distance}, the distance
           along a meridian between the equator and a pole (C{meter}),
           M{b * Elliptic(-e2 / (1 - e2)).E} or M{a * PI / 2}.
        '''
        if self.f:  # complete integral 2nd ...
            # kind: Elliptic(-e2 / (1 - e2)).E
            r = self._elliptic_e22.cE
        else:  # spherical
            r = PI_2
        return Distance(L=self.b * r)

    def Llat(self, lat):
        '''Return the I{meridional length}, the distance along a meridian
           between the equator and a (geodetic) latitude, see C{L}.

           @arg lat: Geodetic latitude (C{degrees90}).

           @return: The meridional length at B{C{lat}}, negative on southern
                    hemisphere (C{meter}).
        '''
        r = self._elliptic_e22.fEd(self.auxParametric(lat)) if self.f else Phi_(lat)
        return Distance(Llat=self.b * r)

    Lmeridian = Llat  # meridional distance

    @Property_RO
    def Mabcd(self):
        '''Get the OSGR meridional coefficients (C{4-Tuple}), C{Airy130} only.
        '''
        if self.n:
            n1, n2, n3 = fpowers(self.n, 3)  # PYCHOK false!
            Mabcd = (fsum_(4, 4 * n1,  5 * n2,  5 * n3) / 4,
                     fsum_(  24 * n1, 24 * n2, 21 * n3) / 8,
                     fsum_(           15 * n2, 15 * n3) / 8,
                                              (35 * n3) / 24)
        else:
            Mabcd = _1_0, _0_0, _0_0, _0_0
        return Mabcd

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{a} or C{Requatorial}.'''
        return self.a

    def m2degrees(self, distance, lat=0):
        '''Convert a distance to an angle along the equator or
           along a parallel of (geodetic) latitude.

           @arg distance: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}) or C{INF} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{distance}} or B{C{lat}}.
       '''
        return degrees(self.m2radians(distance, lat=lat))

    def m2radians(self, distance, lat=0):
        '''Convert a distance to an angle along the equator or
           along a parallel of (geodetic) latitude.

           @arg distance: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{radians}) or C{INF} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{distance}} or B{C{lat}}.
        '''
        r = self.circle4(lat).radius if lat else self.a
        return m2radians(distance, radius=r, lat=0)

    @deprecated_Property_RO
    def minoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{b} or C{Rpolar}.'''
        return self.b

    @Property_RO
    def n(self):
        '''Get the I{3rd flattening} (C{float}), M{f / (2 - f) == (a - b) / (a + b)}, see C{a_b2n}.
        '''
        return self._assert(a_b2n(self.a, self.b), n=f2n(self.f))

    flattening3rd = n

    @deprecated_Property_RO
    def quarteradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{L} or method C{Llat}.'''
        return self.L

    @Property_RO
    def R1(self):
        '''Get the I{mean} earth radius per I{IUGG} (C{meter}), M{(2 * a + b) / 3}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}
                 and method C{Rgeometric}.
        '''
        r = (fsum_(self.a, self.a, self.b) / _3_0) if self.f else self.a
        return Radius(R1=r)

    Rmean = R1

    @Property_RO
    def R2(self):
        '''Get the I{authalic} earth radius (C{meter}), M{sqrt(c2)}.

           @see: C{R2x}, C{c2}, C{area} and U{Earth radius
                 <https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        return Radius(R2=sqrt(self.c2) if self.f else self.a)

    Rauthalic = R2

    @Property_RO
    def R2x(self):
        '''Get the I{authalic} earth radius (C{meter}), M{sqrt(c2x)}.

           @see: C{R2}, C{c2x} and C{areax}.
        '''
        return Radius(R2x=sqrt(self.c2x) if self.f else self.a)

    Rauthalicx = R2x

    @Property_RO
    def R3(self):
        '''Get the I{volumetric} earth radius (C{meter}), M{(a * a * b)**(1/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>} and C{volume}.
        '''
        r = (cbrt(self.b_a) * self.a) if self.f else self.a
        return Radius(R3=r)

    Rvolumetric = R3

    def radians2m(self, rad, lat=0):
        '''Convert an angle to the distance along the equator or
           along a parallel of (geodetic) latitude.

           @arg rad: The angle (C{radians}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Distance (C{meter}, same units as the equatorial
                    and polar radii) or C{0} for near-polar B{C{lat}}.

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{rad}} or B{C{lat}}.
        '''
        r = self.circle4(lat).radius if lat else self.a
        return radians2m(rad, radius=r, lat=0)

    @Property_RO
    def Rbiaxial(self):
        '''Get the I{biaxial, quadratic} mean earth radius (C{meter}), M{sqrt((a**2 + b**2) / 2)}.

           @see: C{Rtriaxial}
        '''
        b = (sqrt((_1_0 + self.b2_a2) * _0_5) * self.a) if self.f else self.a
        return Radius(Rbiaxial=b)

    Requatorial = a  # for consistent naming

    def Rgeocentric(self, lat):
        '''Compute the I{geocentric} earth radius of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Geocentric earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Geocentric Radius
                 <https://WikiPedia.org/wiki/Earth_radius#Geocentric_radius>}
        '''
        r, a = self.a, Phi_(lat)
        if a and self.f:
            if abs(a) < PI_2:
                s2, c2 = _s2_c2(a)
                b2_a2_s2 = self.b2_a2 * s2
                # R == sqrt((a2**2 * c2 + b2**2 * s2) / (a2 * c2 + b2 * s2))
                #   == sqrt(a2**2 * (c2 + (b2 / a2)**2 * s2) / (a2 * (c2 + b2 / a2 * s2)))
                #   == sqrt(a2 * (c2 + (b2 / a2)**2 * s2) / (c2 + (b2 / a2) * s2))
                #   == a * sqrt((c2 + b2_a2 * b2_a2 * s2) / (c2 + b2_a2 * s2))
                #   == a * sqrt((c2 + b2_a2 * b2_a2_s2) / (c2 + b2_a2_s2))
                r *= sqrt((c2 + b2_a2_s2 * self.b2_a2) / (c2 + b2_a2_s2))
            else:
                r  = self.b
        return Radius(Rgeocentric=r)

    @Property_RO
    def Rgeometric(self):
        '''Get the I{geometric} mean earth radius (C{meter}), M{sqrt(a * b)}.

           @see: C{R1}.
        '''
        g = sqrt(self.a * self.b) if self.f else self.a
        return Radius(Rgeometric=g)

    def Rlat(self, lat):
        '''I{Approximate} the earth radius of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Approximate earth radius (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise TypeError: Invalid B{C{lat}}.

           @raise ValueError: Invalid B{C{lat}}.

           @note: C{Rlat(B{90})} equals C{Rpolar}.

           @see: Method C{Rparallel}.
        '''
        # r = a - (a - b) * |lat| / 90
        r = self.a
        if self.f and lat:  # .isEllipsoidal
            r -= (r - self.b) * abs(Lat(lat)) / _90_0
            r  = Radius(Rlat=r)
        return r

    Rpolar = b  # for consistent naming

    @deprecated_Property_RO
    def Rquadratic(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rbiaxial} or C{Rtriaxial}.'''
        return self.Rbiaxial

    @deprecated_Property_RO
    def Rr(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rrectifying}.'''
        return self.Rrectifying

    @Property_RO
    def Rrectifying(self):
        '''Get the I{rectifying} earth radius (C{meter}), M{((a**(3/2) + b**(3/2)) / 2)**(2/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        r = (cbrt2((_1_0 + sqrt3(self.b_a)) * _0_5) * self.a) if self.f else self.a
        return Radius(Rrectifying=r)

    def roc1_(self, sa, ca=None):
        '''Compute the I{prime-vertical}, I{normal} radius of curvature
           of (geodetic) latitude, I{unscaled}.

           @arg sa: Sine of the latitude (C{float}, [-1.0..+1.0]).
           @kwarg ca: Optional cosine of the latitude (C{float}, [-1.0..+1.0])
                      to use an alternate formula.

           @return: The prime-vertical radius of curvature (C{float}).

           @note: The delta between both formulae with C{Ellipsoids.WGS84}
                  is less than 2 nanometer over the entire latitude range.

           @see: Method L{roc2_} and class L{EcefYou}.
        '''
        if not self.f:  # .isSpherical
            n = self.a
        elif ca is None:
            r =  self.e2s2(sa)  # see .roc2_ and _EcefBase._forward
            n = (self.a / sqrt(r)) if r > EPS02 else _0_0
        elif ca:  # derived from EcefYou.forward
            h = hypot(ca, self.b_a * sa) if sa else abs(ca)
            n = self.a / h
        elif sa:
            n = self.a2_b / abs(sa)
        else:
            n = self.a
        return n

    def roc2(self, lat, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal}
           radii of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).
           @kwarg scaled: Scale prime_vertical by C{cos(radians(B{lat}))} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2_} and L{roc1_} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 and meridional and prime vertical U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat), scaled=scaled)

    def roc2_(self, phi, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal}
           radii of curvature of (geodetic) latitude.

           @arg phi: Latitude (C{radians}).
           @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{phi}}.

           @see: Methods L{roc2} and L{roc1_} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 and the meridional and prime vertical U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        a = abs(Phi(phi))
        if self.f:
            r = self.e2s2(sin(a))
            if r > EPS02:
                n = self.a / sqrt(r)
                m = n * self.e12 / r  # PYCHOK attr
            else:
                m = n = _0_0  # PYCHOK attr
        else:
            m = n = self.a
        if scaled and a:
            n *= cos(a) if a < PI_2 else _0_0
        return Curvature2Tuple(Radius(rocMeridional=m),
                               Radius(rocPrimeVertical=n))

    def rocBearing(self, lat, bearing):
        '''Compute the I{directional} radius of curvature
           of (geodetic) latitude and compass direction.

           @arg lat: Latitude (C{degrees90}).
           @arg bearing: Direction (compass C{degrees360}).

           @return: Directional radius of curvature (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{bearing}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        if self.f:
            s2, c2 = _s2_c2(Bearing_(bearing))
            m, n = self.roc2_(Phi_(lat))
            if n < m:  # == n / (c2 * n / m + s2)
                c2 *= n / m
            elif m < n:  # == m / (c2 + s2 * m / n)
                s2 *= m / n
                n = m
            b = n / (c2 + s2)  # == 1 / (c2 / m + s2 / n)
        else:
            b = self.b  # == self.a
        return Radius(rocBearing=b)

    def rocGauss(self, lat):
        '''Compute the I{Gaussian} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Gaussian radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        # using ...
        #    m, n = self.roc2_(Phi_(lat))
        #    return sqrt(m * n)
        # ... requires 1 or 2 sqrt
        if self.f:
            s2, c2 = _s2_c2(Phi_(lat))
            g = self.b / (c2 + self.b2_a2 * s2)
        else:
            g = self.b
        return Radius(rocGauss=g)

    def rocMean(self, lat):
        '''Compute the I{mean} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Mean radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        if self.f:
            m, n = self.roc2_(Phi_(lat))
            m *= 2 * n / (m + n)  # == 2 / (1 / m + 1 / n)
        else:
            m = self.a
        return Radius(rocMean=m)

    def rocMeridional(self, lat):
        '''Compute the I{meridional} radius of curvature of (geodetic) latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Meridional radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2} and L{roc2_} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 and U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).meridional

    rocPolar = a2_b  # synonymous

    def rocPrimeVertical(self, lat):
        '''Compute the I{prime-vertical}, I{normal} radius of curvature
           of (geodetic) latitude, aka the transverse radius of curvature.

           @arg lat: Latitude (C{degrees90}).

           @return: Prime-vertical radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: Methods L{roc2}, L{roc2_} and L{roc1_} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 and U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).prime_vertical

    rocTransverse = rocPrimeVertical  # synonyms

    @deprecated_Property_RO
    def Rs(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Rgeometric}.'''
        return self.Rgeometric

    @Property_RO
    def Rtriaxial(self):
        '''Get the I{triaxial, quadratic} mean earth radius (C{meter}), M{sqrt((3 * a**2 + b**2) / 4)}.

           @see: C{Rbiaxial}
        '''
        t = (sqrt((_3_0 + self.b2_a2) * _0_25) * self.a) if self.f else self.a
        return Radius(Rtriaxial=t)

    def toStr(self, prec=8):  # PYCHOK expected
        '''Return this ellipsoid as a text string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Ellipsoid attributes (C{str}).
        '''
        E = Ellipsoid
        return self._instr(prec, _a_, E.b.name, E.f_.name, _f_, E.f2.name, E.n.name,
                                 E.e.name, E.e2.name, E.e22.name, E.e32.name,
                                 E.A.name, E.L.name, E.R1.name, E.R2.name, E.R3.name,
                                 E.Rbiaxial.name, E.Rtriaxial.name)

    @Property_RO
    def volume(self):
        '''Get the ellipsoid's I{volume} (C{meter**3}), M{4 / 3 * PI * R3**3}.

           @see: C{R3}.
        '''
        return Meter3(volume=self.a2 * self.b * (PI4 / _3_0))


class Ellipsoid2(Ellipsoid):
    '''An L{Ellipsoid} specified by I{equatorial} radius and I{flattening}.
    '''
    def __init__(self, a, f, name=NN):
        '''New L{Ellipsoid2}.

           @arg a: Equatorial radius, semi-axis (C{meter}).
           @arg f: Flattening: (C{float} < 1.0, negative for I{prolate}).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}} or B{C{f}}.

           @note: C{abs(B{f}) < EPS} is forced to C{B{f}=0}, I{spherical}.
                  Negative C{B{f}} produces a I{prolate} ellipsoid.
        '''
        try:
            a = Radius_(a=a, low=EPS,  high=None)  # like Radius_
            f = Float_( f=f, low=None, high=EPS1)
        except (TypeError, ValueError) as x:
            t = instr(self, a=a, f=f, name=name)  # PYCHOK no cover
            raise _ValueError(t, txt=str(x))
        Ellipsoid.__init__(self, a, a_f2b(a, f), name=name)


def _spherical_a_b(a, b):
    '''(INTERNAL) C{True} for spherical or invalid C{a} or C{b}.
    '''
    return a < EPS or b < EPS or abs(a - b) < EPS


def _spherical_f(f):
    '''(INTERNAL) C{True} for spherical or invalid C{f}.
    '''
    return abs(f) < EPS or f > EPS1


def _spherical_f_(f_):
    '''(INTERNAL) C{True} for spherical or invalid C{f_}.
    '''
    return abs(f_) < EPS or abs(f_) > _1_EPS


def a_b2e(a, b):
    '''Return C{e}, the I{1st eccentricity} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The (1st) eccentricity (C{float} or C{0}), M{sqrt(1 - (b / a)**2)}.

       @note: The result is always non-negative and C{0} for I{near-spherical} ellipsoids.
    '''
    return Float(e=sqrt(abs(a_b2e2(a, b))))


def a_b2e2(a, b):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The (1st) eccentricity I{squared} (C{float} or C{0}), M{1 - (b / a)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    return Float(e2=_0_0 if _spherical_a_b(a, b) else (_1_0 - (b / a)**2))


def a_b2e22(a, b):
    '''Return C{e22}, the I{2nd eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The 2nd eccentricity I{squared} (C{float} or C{0}), M{(a / b)**2 - 1}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    return Float(e22=_0_0 if _spherical_a_b(a, b) else ((a / b)**2 - _1_0))


def a_b2e32(a, b):
    '''Return C{e32}, the I{3rd eccentricity squared} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The 3rd eccentricity I{squared} (C{float} or C{0}), M{(a**2 - b**2) / (a**2 + b**2)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    a2, b2 = a**2, b**2
    return Float(e32=_0_0 if _spherical_a_b(a2, b2) else ((a2 - b2) / (a2 + b2)))


def a_b2f(a, b):
    '''Return C{f}, the I{flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The flattening (C{float} or C{0}), M{(a - b) / a}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f = 0 if _spherical_a_b(a, b) else ((a - b) / a)
    return _f_0_0 if _spherical_f(f) else Float(f=f)


def a_b2f_(a, b):
    '''Return C{f_}, the I{inverse flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The inverse flattening (C{float} or C{0}), M{a / (a - b)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    f_ = 0 if _spherical_a_b(a, b) else (a / float(a - b))
    return _f__0_0 if _spherical_f_(f_) else Float(f_=f_)


def a_b2f2(a, b):
    '''Return C{f2}, the I{2nd flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The 2nd flattening (C{float} or C{0}), M{(a - b) / b}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    t = 0 if _spherical_a_b(a, b) else float(a - b)
    return Float(f2=_0_0 if abs(t) < EPS else (t / b))


def a_b2n(a, b):
    '''Return C{n}, the I{3rd flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg b: Polar radius (C{scalar} > 0).

       @return: The 3rd flattening (C{float} or C{0}), M{(a - b) / (a + b)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.
    '''
    t = 0 if _spherical_a_b(a, b) else float(a - b)
    return Float(n=_0_0 if abs(t) < EPS else (t / (a + b)))


def a_f2b(a, f):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{flattening}.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The polar radius (C{float}), M{a * (1 - f)}.
    '''
    b = a if _spherical_f(f) else (a * (_1_0 - f))
    return Radius_(b=a if _spherical_a_b(a, b) else b)


def a_f_2b(a, f_):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{inverse flattening}.

       @arg a: Equatorial radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The polar radius (C{float}), M{a * (f_ - 1) / f_}.
    '''
    b = a if _spherical_f_(f_) else (a * (f_ - _1_0) / f_)
    return Radius_(b=a if _spherical_a_b(a, b) else b)


def b_f2a(b, f):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{flattening}.

       @arg b: Polar radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The equatorial radius (C{float}), M{b / (1 - f)}.
    '''
    t = _1_0 - f
    a = b if abs(t < EPS) else (b / t)
    return Radius_(a=b if _spherical_a_b(a, b) else a)


def b_f_2a(b, f_):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{inverse flattening}.

       @arg b: Polar radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The equatorial radius (C{float}), M{b * f_ / (f_ - 1)}.
    '''
    t = f_ - _1_0
    a = b if _spherical_f_(f_) or abs(t - f_) < EPS \
                               or abs(t) < EPS else (b * f_ / t)
    return Radius_(a=b if _spherical_a_b(a, b) else a)


def f2e2(f):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The (1st) eccentricity I{squared} (C{float} < 1), M{f * (2 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(e2=_0_0 if _spherical_f(f) else (f * (_2_0 - f)))


def f2e22(f):
    '''Return C{e22}, the I{2nd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The 2nd eccentricity I{squared} (C{float} > -1 or C{INF}),
                M{f * (2 - f) / (1 - f)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for near-spherical ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (1 - e2) == f * (2 - f) / (1 - f)**2
    t = (_1_0 - f)**2
    return Float(e22=INF if t < EPS else (f2e2(f) / t))  # PYCHOK type


def f2e32(f):
    '''Return C{e32}, the I{3rd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The 3rd eccentricity I{squared} (C{float}), M{f * (2 - f) / (1 + (1 - f)**2)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (2 - e2) == f * (2 - f) / (1 + (1 - f)**2)
    e2 = f2e2(f)
    return Float(e32=e2 / (_2_0 - e2))


def f_2f(f_):
    '''Return C{f}, the I{flattening} for a given I{inverse flattening}.

       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The flattening (C{float} or C{0}), M{1 / f_}.

       @note: The result is positive for I{oblate}, negative for I{prolate}
              or C{0} for I{near-spherical} ellipsoids.
    '''
    f = 0 if _spherical_f_(f_) else _1_0 / f_
    return _f_0_0 if _spherical_f(f) else Float(f=f)  # PYCHOK type


def f2f_(f):
    '''Return C{f_}, the I{inverse flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The inverse flattening (C{float} or C{0}), M{1 / f}.

       @note: The result is positive for I{oblate}, negative for I{prolate}
              or C{0} for I{near-spherical} ellipsoids.
    '''
    f_ = 0 if _spherical_f(f) else _1_0 / f
    return _f__0_0 if _spherical_f_(f_) else Float(f_=f_)  # PYCHOK type


def f2f2(f):
    '''Return C{f2}, the I{2nd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The 2nd flattening (C{float} or C{INF}), M{f / (1 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate}
              or C{0} for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = _1_0 - f
    return Float(f2=_0_0 if _spherical_f(f) else
                    (INF if  abs(t) < EPS else (f / t)))  # PYCHOK type


def f2n(f):
    '''Return C{n}, the I{3rd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1, negative for I{prolate}).

       @return: The 3rd flattening (-1 <= C{float} < 1), M{f / (2 - f)}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(n=_0_0 if _spherical_f(f) else (f / float(_2_0 - f)))


def n2e2(n):
    '''Return C{e2}, the I{1st eccentricity squared} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The (1st) eccentricity I{squared} (C{float} or -INF), M{4 * n / (1 + n)**2}.

       @note: The result is positive for I{oblate}, negative for I{prolate} or C{0}
              for I{near-spherical} ellipsoids.

       @see: U{Flattening<https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = (_1_0 + n)**2
    return Float(e2=_0_0 if abs(n) < EPS else
                   (-INF if     t  < EPS else (_4_0 * n / t)))


def _normal2(px, py, E):
    '''(INTERNAL) Nearest point on a 2-D ellipse.
    '''
    a, b = E.a, E.b
    if min(px, py, a, b) < EPS0:
        raise _AssertionError(px=px, py=py, a=a, b=b, E=E)

    a2 = a - b * E.b_a
    b2 = b - a * E.a_b
    tx = ty = 0.70710678118
    for _ in range(9):  # 4..5 trips max
        ex = a2 * tx**3
        ey = b2 * ty**3

        qx = px - ex
        qy = py - ey
        q  = hypot(qx, qy)
        if q < EPS0:
            break
        r = hypot(ex - tx * a, ey - ty * b) / q

        sx, tx = tx, min(_1_0, max(0, (ex + qx * r) / a))
        sy, ty = ty, min(_1_0, max(0, (ey + qy * r) / b))
        t = hypot(ty, tx)
        if t < EPS0:
            break
        tx = tx / t  # /= t chokes PyChecker
        ty = ty / t
        if max(abs(sx - tx), abs(sy - ty)) < EPS:
            break

    tx *= a / px
    ty *= b / py
    return tx, ty  # x and y as fractions


def n2f(n):
    '''Return C{f}, the I{flattening} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The flattening (C{float} or -INF), M{2 * n / (1 + n)}.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = n + _1_0
    f = 0 if abs(n) < EPS else (-INF if t < EPS else _2_0 * n / t)
    return _f_0_0 if _spherical_f(f) else Float(f=f)


def n2f_(n):
    '''Return C{f_}, the I{inverse flattening} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} < 1).

       @return: The inverse flattening (C{float} or C{0}), M{1 / f}.

       @see: L{n2f} and L{f2f_}.
    '''
    return f2f_(n2f(n))


class Ellipsoids(_NamedEnum):
    '''(INTERNAL) L{Ellipsoid} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, a, b, f_, **kwds):
        '''(INTERNAL) Instantiate the L{Ellipsoid}.
        '''
        return Ellipsoid(a, b, f_, **kwds)

Ellipsoids = Ellipsoids(Ellipsoid)  # PYCHOK singleton
'''Some pre-defined L{Ellipsoid}s, all I{lazily} instantiated.'''
# <https://www.GNU.org/software/gama/manual/html_node/Supported-ellipsoids.html>
# <https://w3.Energistics.org/archive/Epicentre/Epicentre_v3.0/DataModel/
#         LogicalDictionary/StandardValues/ellipsoid.html>
# <https://kb.OSU.edu/dspace/handle/1811/77986>
Ellipsoids._assert(  # <https://WikiPedia.org/wiki/Earth_ellipsoid>
    Airy1830       = _lazy(_Airy1830_,       *_T(6377563.396, _0_0,               299.3249646)),   # b=6356256.909
    AiryModified   = _lazy(_AiryModified_,   *_T(6377340.189, _0_0,               299.3249646)),  # b=6356034.448
#   ANS            = _lazy('ANS',            *_T(6378160.0,   _0_0,               298.25)),  # Australian Nat. Spheroid
    Australia1966  = _lazy('Australia1966',  *_T(6378160.0,   _0_0,               298.25)),  # b=6356774.7192
#   Bessel1841     = _lazy(_Bessel1841_,     *_T(6377397.155,  6356078.963,       299.152815351)),
    Bessel1841     = _lazy(_Bessel1841_,     *_T(6377397.155,  6356078.962818,    299.152812797)),
    Clarke1866     = _lazy(_Clarke1866_,     *_T(6378206.4,    6356583.8,         294.978698214)),
    Clarke1880     = _lazy('Clarke1880',     *_T(6378249.145,  6356514.86954978,  293.465)),
    Clarke1880IGN  = _lazy(_Clarke1880IGN_,  *_T(6378249.2,    6356515.0,         293.466021294)),
    Clarke1880Mod  = _lazy('Clarke1880Mod',  *_T(6378249.145,  6356514.96582849,  293.4663)),
    CPM1799        = _lazy('CPM1799',        *_T(6375738.7,    6356671.92557493,  334.39)),  # Comm. des Poids et Mesures
    Delambre1810   = _lazy('Delambre1810',   *_T(6376428.0,    6355957.92616372,  311.5)),  # Belgium
    Engelis1985    = _lazy('Engelis1985',    *_T(6378136.05,   6356751.32272154,  298.2566)),
    Everest1969    = _lazy('Everest1969',    *_T(6377295.664,  6356094.667915,    300.801699997)),
    Fisher1968     = _lazy('Fisher1968',     *_T(6378150.0,    6356768.33724438,  298.3)),
    GEM10C         = _lazy('GEM10C',         *_T(R_MA,         6356752.31424783,  298.2572236)),
#   GPES           = _lazy('GPES',           *_T(6378135.0,    6356750.0,        _0_0)),  # "Gen. Purpose Earth Spheroid"
    GRS67          = _lazy('GRS67',          *_T(6378160.0,   _0_0,               298.247167427)),  # Lucerne b=6356774.516
    GRS80          = _lazy(_GRS80_,          *_T(R_MA,         6356752.314140347, 298.257222101)),  # ITRS, ETRS89
#   Hayford1924    = _lazy('Hayford1924',    *_T(6378388.0,    6356911.94612795, _0_0)),  # aka Intl1924 f_=297
    Helmert1906    = _lazy('Helmert1906',    *_T(6378200.0,    6356818.16962789,  298.3)),
    IERS1989       = _lazy('IERS1989',       *_T(6378136.0,   _0_0,               298.257)),  # b=6356751.302
    IERS1992TOPEX  = _lazy('IERS1992TOPEX',  *_T(6378136.3,    6356751.61659215,  298.257223563)),  # IERS/TOPEX/Poseidon/McCarthy
    IERS2003       = _lazy('IERS2003',       *_T(6378136.6,    6356751.85797165,  298.25642)),
    Intl1924       = _lazy(_Intl1924_,       *_T(6378388.0,   _0_0,               297.0)),  # aka Hayford b=6356911.9462795
    Intl1967       = _lazy('Intl1967',       *_T(6378157.5,    6356772.2,         298.24961539)),  # New Int'l
    Krassovski1940 = _lazy(_Krassovski1940_, *_T(6378245.0,    6356863.01877305,  298.3)),  # spelling
    Krassowsky1940 = _lazy(_Krassowsky1940_, *_T(6378245.0,    6356863.01877305,  298.3)),  # spelling
    Maupertuis1738 = _lazy('Maupertuis1738', *_T(6397300.0,    6363806.28272251,  191.0)),  # France
    Mercury1960    = _lazy('Mercury1960',    *_T(6378166.0,    6356784.28360711,  298.3)),
    Mercury1968Mod = _lazy('Mercury1968Mod', *_T(6378150.0,    6356768.33724438,  298.3)),
    NWL1965        = _lazy('NWL1965',        *_T(6378145.0,    6356759.76948868,  298.25)),  # Naval Weapons Lab.
    OSU86F         = _lazy('OSU86F',         *_T(6378136.2,    6356751.51693008,  298.2572236)),
    OSU91A         = _lazy('OSU91A',         *_T(6378136.3,    6356751.6165948,   298.2572236)),
#   Plessis1817    = _lazy('Plessis1817',    *_T(6397523.0,    6355863.0,         153.56512242)),  # XXX incorrect?
    Plessis1817    = _lazy('Plessis1817',    *_T(6376523.0,    6355862.93325557,  308.64)),  # XXX IGN France 1972
    SGS85          = _lazy('SGS85',          *_T(6378136.0,    6356751.30156878,  298.257)),  # Soviet Geodetic System
    SoAmerican1969 = _lazy('SoAmerican1969', *_T(6378160.0,    6356774.71919531,  298.25)),  # South American
    Struve1860     = _lazy('Struve1860',     *_T(6378298.3,    6356657.14266956,  294.73)),
    WGS60          = _lazy('WGS60',          *_T(6378165.0,    6356783.28695944,  298.3)),
    WGS66          = _lazy('WGS66',          *_T(6378145.0,    6356759.76948868,  298.25)),
    WGS72          = _lazy(_WGS72_,          *_T(6378135.0,   _0_0,               298.26)),  # b=6356750.52
    WGS84          = _lazy(_WGS84_,          *_T(R_MA,        _0_0,              _f_WGS84)),  # GPS b=6356752.3142451793
#   Prolate        = _lazy('Prolate',        *_T(6356752.3,    R_MA,             _0_0)),
    Sphere         = _lazy(_Sphere_,         *_T(R_M,          R_M,              _0_0)),  # pseudo
    SphereAuthalic = _lazy('SphereAuthalic', *_T(R_FM,         R_FM,             _0_0)),  # pseudo
    SpherePopular  = _lazy('SpherePopular',  *_T(R_MA,         R_MA,             _0_0))   # EPSG:3857 Spheroid
)


if __name__ == '__main__':

    from pygeodesy.interns import _COMMA_, _NL_, _NL_hash_, _NL_var_
    from pygeodesy.named import nameof

    for E in (Ellipsoids.WGS84, Ellipsoids.GRS80,  # NAD83,
              Ellipsoids.Sphere, Ellipsoids.SpherePopular,
              Ellipsoid(Ellipsoids.WGS84.b, Ellipsoids.WGS84.a, name='_Prolate')):
        e = f2n(E.f) - E.n
        t = NN(_COMMA_, _NL_hash_, _SPACE_)(E.toStr(prec=10),  # re-callable
               'e=%s, f_=%s, f=%s, n=%s (%s)' % (fstr(E.e,  prec=13, fmt=Fmt.e),
                                                 fstr(E.f_, prec=13, fmt=Fmt.e),
                                                 fstr(E.f,  prec=13, fmt=Fmt.e),
                                                 fstr(E.n,  prec=13, fmt=Fmt.e),
                                                 fstr(e,    prec=9,  fmt=Fmt.e),),
               '%s=(%s)'   % (Ellipsoid.AlphaKs.name, fstr(E.AlphaKs, prec=20),),
               '%s= (%s)'  % (Ellipsoid.BetaKs.name,  fstr(E.BetaKs,  prec=20),),
               '%s= %s'    % (nameof(Ellipsoid.KsOrder),   E.KsOrder),  # property
               '%s=  (%s)' % (Ellipsoid.Mabcd.name,   fstr(E.Mabcd,   prec=20),))
        print('%s%s: %s' % (_NL_hash_, _DOT_(E.classname, E.name), t))

    # __doc__ of this file, force all into registry
    t = [NN] + Ellipsoids.toRepr(all=True).split(_NL_)
    print(_NL_var_.join(i.strip(_COMMA_) for i in t))

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

# % python3 -m pygeodesy.ellipsoids

# Ellipsoid.WGS84: name='WGS84', a=6378137, b=6356752.3142451793, f_=298.257223563, f=0.0033528107, f2=0.0033640898, n=0.0016792204, e=0.0818191908, e2=0.00669438, e22=0.0067394967, e32=0.0033584313, A=6367449.1458234144, L=10001965.7293127235, R1=6371008.7714150595, R2=6371007.1809184747, R3=6371000.7900091587, Rbiaxial=6367453.6345163295, Rtriaxial=6372797.5559594007,
#  e=8.1819190842622e-02, f_=2.98257223563e+02, f=3.3528106647475e-03, n=1.6792203863837e-03 (0.0e+00),
#  AlphaKs=(0.00083773182062447786, 0.00000076085277735726, 0.00000000119764550324, 0.00000000000242917068, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#  BetaKs= (0.00083773216405795667, 0.0000000590587015222, 0.00000000016734826653, 0.00000000000021647981, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#  KsOrder= 8,
#  Mabcd=  (1.00168275103155868244, 0.00504613293193333871, 0.00000529596776243457, 0.00000000690525779769)

# Ellipsoid.GRS80: name='GRS80', a=6378137, b=6356752.3141403468, f_=298.257222101, f=0.0033528107, f2=0.0033640898, n=0.0016792204, e=0.081819191, e2=0.00669438, e22=0.0067394968, e32=0.0033584313, A=6367449.1457710434, L=10001965.7292304579, R1=6371008.7713801153, R2=6371007.1808835147, R3=6371000.7899741363, Rbiaxial=6367453.6344640013, Rtriaxial=6372797.5559332585,
#  e=8.1819191042833e-02, f_=2.98257222101e+02, f=3.3528106811837e-03, n=1.6792203946295e-03 (0.0e+00),
#  AlphaKs=(0.00083773182472890429, 0.00000076085278481561, 0.00000000119764552086, 0.00000000000242917073, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#  BetaKs= (0.0008377321681623882, 0.00000005905870210374, 0.000000000167348269, 0.00000000000021647982, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#  KsOrder= 8,
#  Mabcd=  (1.00168275103983916985, 0.0050461329567537995, 0.00000529596781448937, 0.00000000690525789941)

# Ellipsoid.Sphere: name='Sphere', a=6371008.7714149999, b=6371008.7714149999, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, A=6371008.7714149999, L=10007557.1761167478, R1=6371008.7714149999, R2=6371008.7714149999, R3=6371008.7714149999, Rbiaxial=6371008.7714149999, Rtriaxial=6371008.7714149999,
#  e=0.0e+00, f_=0.0e+00, f=0.0e+00, n=0.0e+00 (0.0e+00),
#  AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#  BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#  KsOrder= 8,
#  Mabcd=  (1.0, 0.0, 0.0, 0.0)

# Ellipsoid.SpherePopular: name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, A=6378137, L=10018754.171394622, R1=6378137, R2=6378137, R3=6378137, Rbiaxial=6378137, Rtriaxial=6378137,
#  e=0.0e+00, f_=0.0e+00, f=0.0e+00, n=0.0e+00 (0.0e+00),
#  AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#  BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#  KsOrder= 8,
#  Mabcd=  (1.0, 0.0, 0.0, 0.0)

# Ellipsoid._Prolate: name='_Prolate', a=6356752.3142451793, b=6378137, f_=-297.257223563, f=-0.0033640898, f2=-0.0033528107, n=-0.0016792204, e=0.0820944379, e2=-0.0067394967, e22=-0.00669438, e32=-0.0033584313, A=6367449.1458234154, L=10035500.5204500332, R1=6363880.5428301198, R2=6363878.9413582645, R3=6363872.5644020075, Rbiaxial=6367453.6345163304, Rtriaxial=6362105.2243882548,
#  e=8.2094437949696e-02, f_=-2.97257223563e+02, f=-3.3640898209765e-03, n=-1.6792203863837e-03 (0.0e+00),
#  AlphaKs=(-0.00084149152514366627, 0.00000076653480614871, -0.00000000120934503389, 0.0000000000024576225, -0.00000000000000578863, 0.00000000000000001502, -0.00000000000000000004, 0.0),
#  BetaKs= (-0.00084149187224351817, 0.00000005842735196773, -0.0000000001680487236, 0.00000000000021706261, -0.00000000000000038002, 0.00000000000000000073, -0.0, 0.0),
#  KsOrder= 8,
#  Mabcd=  (0.99832429842120640195, -0.00502921424529705757, 0.00000527821138524052, -0.00000000690525779769)
