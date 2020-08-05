
# -*- coding: utf-8 -*-

u'''Classes L{Datum}, L{Ellipsoid} and L{Transform} and registries thereof
L{Datums}, L{Ellipsoids} and L{Transforms}, respectively.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
including datums and ellipsoid parameters for different geographic coordinate
systems and methods for converting between them and to cartesian coordinates.
Transcribed from JavaScript originals by I{(C) Chris Veness 2005-2016} and
published under the same MIT Licence**, see U{latlon-ellipsoidal.js
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

Historical geodetic datums: a latitude/longitude point defines a geographic
location on or above/below the earth’s surface, measured in degrees from
the equator and the International Reference Meridian and meters above the
ellipsoid, and based on a given datum.  The datum is based on a reference
ellipsoid and tied to geodetic survey reference points.

Modern geodesy is generally based on the WGS84 datum (as used for instance
by GPS systems), but previously various reference ellipsoids and datum
references were used.

The UK Ordnance Survey National Grid References are still based on the otherwise
historical OSGB36 datum, q.v. U{Ordnance Survey 'A guide to coordinate systems
in Great Britain', Section 6<https://www.OrdnanceSurvey.co.UK/docs/support/
guide-coordinate-systems-great-britain.pdf>}.

@newfield example: Example, Examples

@var Datums.BD72: Datum(name='BD72', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.BD72)
@var Datums.DHDN: Datum(name='DHDN', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.DHDN)
@var Datums.ED50: Datum(name='ED50', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
@var Datums.GDA2020: Datum(name='GDA2020', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.GRS80: Datum(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.Irl1975: Datum(name='Irl1975', ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
@var Datums.Krassovski1940: Datum(name='Krassovski1940', ellipsoid=Ellipsoids.Krassovski1940, transform=Transforms.Krassovski1940)
@var Datums.Krassowsky1940: Datum(name='Krassowsky1940', ellipsoid=Ellipsoids.Krassowsky1940, transform=Transforms.Krassowsky1940)
@var Datums.MGI: Datum(name='MGI', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.MGI)
@var Datums.NAD27: Datum(name='NAD27', ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
@var Datums.NAD83: Datum(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
@var Datums.NTF: Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
@var Datums.OSGB36: Datum(name='OSGB36', ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
@var Datums.Potsdam: Datum(name='Potsdam', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.Bessel1841)
@var Datums.Sphere: Datum(name='Sphere', ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
@var Datums.TokyoJapan: Datum(name='TokyoJapan', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
@var Datums.WGS72: Datum(name='WGS72', ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
@var Datums.WGS84: Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)

@var Ellipsoids.Airy1830: Ellipsoid(name='Airy1830', a=6377563.396, b=6356256.909, f_=299.3249646, f=0.003340851, e=0.081673374, e2=0.00667054, e12=0.99332946, e22=0.006715335, n=0.00167322, R1=6370461.233666666, R2=6370459.654589443, R3=6370453.309866445, Rr=6366914.608805892, Rs=6366901.239881964)
@var Ellipsoids.AiryModified: Ellipsoid(name='AiryModified', a=6377340.189, b=6356034.448, f_=299.3249646, f=0.003340851, e=0.081673374, e2=0.00667054, e12=0.99332946, e22=0.006715335, n=0.00167322, R1=6370238.275333334, R2=6370236.696361165, R3=6370230.351810658, Rr=6366691.774649803, Rs=6366678.406194146)
@var Ellipsoids.Australia1966: Ellipsoid(name='Australia1966', a=6378160, b=6356774.719, f_=298.25, f=0.003352892, e=0.08182018, e2=0.006694542, e12=0.993305458, e22=0.006739661, n=0.001679261, R1=6371031.573, R2=6371029.982388151, R3=6371023.591178183, Rr=6367471.848433915, Rs=6367458.38162583)
@var Ellipsoids.Bessel1841: Ellipsoid(name='Bessel1841', a=6377397.155, b=6356078.962818, f_=299.1528128, f=0.003342773, e=0.081696831, e2=0.006674372, e12=0.993325628, e22=0.006719219, n=0.001674185, R1=6370291.090939333, R2=6370289.510126558, R3=6370283.158215224, Rr=6366742.520233163, Rs=6366729.136254413)
@var Ellipsoids.CPM1799: Ellipsoid(name='CPM1799', a=6375738.7, b=6356671.92557493, f_=334.39, f=0.00299052, e=0.077279343, e2=0.005972097, e12=0.994027903, e22=0.006007977, n=0.001497499, R1=6369383.108524977, R2=6369381.843415799, R3=6369376.762470212, Rr=6366208.881847335, Rs=6366198.174663714)
@var Ellipsoids.Clarke1866: Ellipsoid(name='Clarke1866', a=6378206.4, b=6356583.8, f_=294.978698214, f=0.003390075, e=0.082271854, e2=0.006768658, e12=0.993231342, e22=0.006814785, n=0.001697916, R1=6370998.866666667, R2=6370997.240632997, R3=6370990.706598808, Rr=6367399.689168951, Rs=6367385.921655473)
@var Ellipsoids.Clarke1880: Ellipsoid(name='Clarke1880', a=6378249.145, b=6356514.86954978, f_=293.465, f=0.003407561, e=0.0824834, e2=0.006803511, e12=0.993196489, e22=0.006850116, n=0.001706689, R1=6371004.386516593, R2=6371002.743669633, R3=6370996.141916499, Rr=6367386.643979662, Rs=6367372.733858579)
@var Ellipsoids.Clarke1880IGN: Ellipsoid(name='Clarke1880IGN', a=6378249.2, b=6356515, f_=293.466021294, f=0.00340755, e=0.082483257, e2=0.006803488, e12=0.993196512, e22=0.006850092, n=0.001706683, R1=6371004.466666666, R2=6371002.823831111, R3=6370996.22212394, Rr=6367386.736672514, Rs=6367372.826648208)
@var Ellipsoids.Clarke1880Mod: Ellipsoid(name='Clarke1880Mod', a=6378249.145, b=6356514.96582849, f_=293.4663, f=0.003407546, e=0.082483218, e2=0.006803481, e12=0.993196519, e22=0.006850086, n=0.001706681, R1=6371004.418609496, R2=6371002.775777078, R3=6370996.174082516, Rr=6367386.692077903, Rs=6367372.782080163)
@var Ellipsoids.Delambre1810: Ellipsoid(name='Delambre1810', a=6376428, b=6355957.92616372, f_=311.5, f=0.003210273, e=0.080063974, e2=0.00641024, e12=0.99358976, e22=0.006451596, n=0.001607717, R1=6369604.642054573, R2=6369603.184197493, R3=6369597.327390675, Rr=6366197.076842674, Rs=6366184.735554905)
@var Ellipsoids.Engelis1985: Ellipsoid(name='Engelis1985', a=6378136.05, b=6356751.32272154, f_=298.2566, f=0.003352818, e=0.081819276, e2=0.006694394, e12=0.993305606, e22=0.006739511, n=0.001679224, R1=6371007.807573847, R2=6371006.217070852, R3=6370999.826135725, Rr=6367448.175078915, Rs=6367434.70891814)
@var Ellipsoids.Everest1969: Ellipsoid(name='Everest1969', a=6377295.664, b=6356094.667915, f_=300.8017, f=0.003324449, e=0.081472981, e2=0.006637847, e12=0.993362153, e22=0.006682202, n=0.001664992, R1=6370228.665305, R2=6370227.10178534, R3=6370220.819516171, Rr=6366699.578394239, Rs=6366686.341077896)
@var Ellipsoids.Fisher1968: Ellipsoid(name='Fisher1968', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371022.77908146, R2=6371021.189037352, R3=6371014.799950343, Rr=6367463.656043013, Rs=6367450.193774211)
@var Ellipsoids.GEM10C: Ellipsoid(name='GEM10C', a=6378137, b=6356752.31424783, f_=298.2572236, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371008.771415944, R2=6371007.180919358, R3=6371000.790010039, Rr=6367449.145823942, Rs=6367435.679717519)
@var Ellipsoids.GRS67: Ellipsoid(name='GRS67', a=6378160, b=6356774.516, f_=298.247167427, f=0.003352924, e=0.081820568, e2=0.006694605, e12=0.993305395, e22=0.006739725, n=0.001679277, R1=6371031.505333333, R2=6371029.914708731, R3=6371023.523359839, Rr=6367471.74701921, Rs=6367458.279955242)
@var Ellipsoids.GRS80: Ellipsoid(name='GRS80', a=6378137, b=6356752.314140347, f_=298.257222101, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371008.771380115, R2=6371007.180883513, R3=6371000.789974131, Rr=6367449.145770245, Rs=6367435.679663688)
@var Ellipsoids.Helmert1906: Ellipsoid(name='Helmert1906', a=6378200, b=6356818.16962789, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371072.723209296, R2=6371071.133152724, R3=6371064.744015628, Rr=6367513.572269944, Rs=6367500.10989561)
@var Ellipsoids.IERS1989: Ellipsoid(name='IERS1989', a=6378136, b=6356751.302, f_=298.257, f=0.003352813, e=0.081819221, e2=0.006694385, e12=0.993305615, e22=0.006739502, n=0.001679222, R1=6371007.767333333, R2=6371006.176906484, R3=6370999.785917024, Rr=6367448.139705879, Rs=6367434.673581903)
@var Ellipsoids.IERS1992TOPEX: Ellipsoid(name='IERS1992TOPEX', a=6378136.3, b=6356751.61659215, f_=298.257223563, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371008.072197382, R2=6371006.481700972, R3=6371000.090792353, Rr=6367448.44699561, Rs=6367434.980890662)
@var Ellipsoids.IERS2003: Ellipsoid(name='IERS2003', a=6378136.6, b=6356751.85797165, f_=298.25642, f=0.00335282, e=0.081819301, e2=0.006694398, e12=0.993305602, e22=0.006739515, n=0.001679225, R1=6371008.352657217, R2=6371006.762152168, R3=6371000.371208764, Rr=6367448.717709779, Rs=6367435.251531576)
@var Ellipsoids.Intl1924: Ellipsoid(name='Intl1924', a=6378388, b=6356911.946, f_=297, f=0.003367003, e=0.08199189, e2=0.00672267, e12=0.99327733, e22=0.00676817, n=0.001686341, R1=6371229.315333334, R2=6371227.711270464, R3=6371221.265832124, Rr=6367654.499992854, Rs=6367640.919007843)
@var Ellipsoids.Intl1967: Ellipsoid(name='Intl1967', a=6378157.5, b=6356772.2, f_=298.24961539, f=0.003352896, e=0.081820233, e2=0.00669455, e12=0.99330545, e22=0.00673967, n=0.001679263, R1=6371029.066666666, R2=6371027.476083895, R3=6371021.084827519, Rr=6367469.338943659, Rs=6367455.872106339)
@var Ellipsoids.Krassovski1940: Ellipsoid(name='Krassovski1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.693674386, Rr=6367558.496874185, Rs=6367545.034404869)
@var Ellipsoids.Krassowsky1940: Ellipsoid(name='Krassowsky1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.693674386, Rr=6367558.496874185, Rs=6367545.034404869)
@var Ellipsoids.Maupertuis1738: Ellipsoid(name='Maupertuis1738', a=6397300, b=6363806.28272251, f_=191, f=0.005235602, e=0.102194876, e2=0.010443793, e12=0.989556207, e22=0.010554017, n=0.002624672, R1=6386135.42757417, R2=6386131.541448465, R3=6386115.886282292, Rr=6380564.130113637, Rs=6380531.163818629)
@var Ellipsoids.Mercury1960: Ellipsoid(name='Mercury1960', a=6378166, b=6356784.28360711, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371038.76120237, R2=6371037.171154276, R3=6371030.782051236, Rr=6367479.629235635, Rs=6367466.166933062)
@var Ellipsoids.Mercury1968Mod: Ellipsoid(name='Mercury1968Mod', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371022.77908146, R2=6371021.189037352, R3=6371014.799950343, Rr=6367463.656043013, Rs=6367450.193774211)
@var Ellipsoids.NWL1965: Ellipsoid(name='NWL1965', a=6378145, b=6356759.76948868, f_=298.25, f=0.003352892, e=0.08182018, e2=0.006694542, e12=0.993305458, e22=0.006739661, n=0.001679261, R1=6371016.58982956, R2=6371014.999254003, R3=6371008.60802666, Rr=6367456.873667616, Rs=6367443.406891448)
@var Ellipsoids.OSU86F: Ellipsoid(name='OSU86F', a=6378136.2, b=6356751.51693008, f_=298.2572236, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371007.972310026, R2=6371006.381813641, R3=6370999.990905124, Rr=6367448.347164505, Rs=6367434.88105977)
@var Ellipsoids.OSU91A: Ellipsoid(name='OSU91A', a=6378136.3, b=6356751.6165948, f_=298.2572236, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371008.072198267, R2=6371006.481701856, R3=6371000.090793238, Rr=6367448.446996935, Rs=6367434.98089199)
@var Ellipsoids.Plessis1817: Ellipsoid(name='Plessis1817', a=6376523, b=6355862.93325557, f_=308.64, f=0.003240021, e=0.080433474, e2=0.006469544, e12=0.993530456, e22=0.006511671, n=0.001622639, R1=6369636.31108519, R2=6369634.826085826, R3=6369628.859996674, Rr=6366197.15710669, Rs=6366184.585664447)
@var Ellipsoids.SGS85: Ellipsoid(name='SGS85', a=6378136, b=6356751.30156878, f_=298.257, f=0.003352813, e=0.081819221, e2=0.006694385, e12=0.993305615, e22=0.006739502, n=0.001679222, R1=6371007.767189593, R2=6371006.176690875, R3=6370999.785772962, Rr=6367448.13949045, Rs=6367434.673365931)
@var Ellipsoids.SoAmerican1969: Ellipsoid(name='SoAmerican1969', a=6378160, b=6356774.71919531, f_=298.25, f=0.003352892, e=0.08182018, e2=0.006694542, e12=0.993305458, e22=0.006739661, n=0.001679261, R1=6371031.573065103, R2=6371029.982485807, R3=6371023.591243432, Rr=6367471.848531488, Rs=6367458.38172365)
@var Ellipsoids.Sphere: Ellipsoid(name='Sphere', a=6371008.771415, b=6371008.771415, f_=0, f=0, e=0, e2=0, e12=1, e22=0, n=0, R1=6371008.771415, R2=6371008.771415, R3=6371008.771415, Rr=6371008.771415, Rs=6371008.771415)
@var Ellipsoids.SphereAuthalic: Ellipsoid(name='SphereAuthalic', a=6371000, b=6371000, f_=0, f=0, e=0, e2=0, e12=1, e22=0, n=0, R1=6371000, R2=6371000, R3=6371000, Rr=6371000, Rs=6371000)
@var Ellipsoids.SpherePopular: Ellipsoid(name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, e=0, e2=0, e12=1, e22=0, n=0, R1=6378137, R2=6378137, R3=6378137, Rr=6378137, Rs=6378137)
@var Ellipsoids.Struve1860: Ellipsoid(name='Struve1860', a=6378298.3, b=6356657.14266956, f_=294.73, f=0.003392936, e=0.082306499, e2=0.00677436, e12=0.99322564, e22=0.006820565, n=0.001699351, R1=6371084.580889854, R2=6371082.952089875, R3=6371076.406914177, Rr=6367482.318324656, Rs=6367468.527348378)
@var Ellipsoids.WGS60: Ellipsoid(name='WGS60', a=6378165, b=6356783.28695944, f_=298.3, f=0.00335233, e=0.081813334, e2=0.006693422, e12=0.993306578, e22=0.006738525, n=0.001678979, R1=6371037.762319813, R2=6371036.172271968, R3=6371029.783169931, Rr=6367478.630911094, Rs=6367465.168610634)
@var Ellipsoids.WGS66: Ellipsoid(name='WGS66', a=6378145, b=6356759.76948868, f_=298.25, f=0.003352892, e=0.08182018, e2=0.006694542, e12=0.993305458, e22=0.006739661, n=0.001679261, R1=6371016.58982956, R2=6371014.999254003, R3=6371008.60802666, Rr=6367456.873667616, Rs=6367443.406891448)
@var Ellipsoids.WGS72: Ellipsoid(name='WGS72', a=6378135, b=6356750.52, f_=298.26, f=0.003352779, e=0.081818811, e2=0.006694318, e12=0.993305682, e22=0.006739434, n=0.001679205, R1=6371006.84, R2=6371005.249530816, R3=6370998.858745317, Rr=6367447.248614991, Rs=6367433.78276368)
@var Ellipsoids.WGS84: Ellipsoid(name='WGS84', a=6378137, b=6356752.31425, f_=298.257223563, f=0.003352811, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.006739497, n=0.00167922, R1=6371008.771416667, R2=6371007.180920884, R3=6371000.790010764, Rr=6367449.145825027, Rs=6367435.679718607)

@var Transforms.BD72: Transform(name='BD72', tx=106.86863, ty=-52.29778, tz=103.72389, rx=-0, ry=-0, rz=-0.00001, s=1.2727, s1=1, sx=-0.33657, sy=-0.45696, sz=-1.84218)
@var Transforms.Bessel1841: Transform(name='Bessel1841', tx=-582, ty=-105, tz=-414, rx=-0.00001, ry=-0, rz=0.00001, s=-8.3, s1=0.99999, sx=-1.04, sy=-0.35, sz=3.08)
@var Transforms.Clarke1866: Transform(name='Clarke1866', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.DHDN: Transform(name='DHDN', tx=-591.28, ty=-81.35, tz=-396.39, rx=0.00001, ry=-0, rz=-0.00001, s=-9.82, s1=0.99999, sx=1.477, sy=-0.0736, sz=-1.458)
@var Transforms.ED50: Transform(name='ED50', tx=89.5, ty=93.8, tz=123.1, rx=0, ry=0, rz=0, s=-1.2, s1=1, sx=0, sy=0, sz=0.156)
@var Transforms.Irl1965: Transform(name='Irl1965', tx=-482.53, ty=130.596, tz=-564.557, rx=0.00001, ry=0, rz=0, s=-8.15, s1=0.99999, sx=1.042, sy=0.214, sz=0.631)
@var Transforms.Irl1975: Transform(name='Irl1975', tx=-482.53, ty=130.596, tz=-564.557, rx=-0.00001, ry=-0, rz=-0, s=-1.1, s1=1, sx=-1.042, sy=-0.214, sz=-0.631)
@var Transforms.Krassovski1940: Transform(name='Krassovski1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.Krassowsky1940: Transform(name='Krassowsky1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.MGI: Transform(name='MGI', tx=-577.326, ty=-90.129, tz=-463.92, rx=0.00002, ry=0.00001, rz=0.00003, s=-2.423, s1=1, sx=5.137, sy=1.474, sz=5.297)
@var Transforms.NAD27: Transform(name='NAD27', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.NAD83: Transform(name='NAD83', tx=1.004, ty=-1.91, tz=-0.515, rx=0, ry=0, rz=0, s=-0.0015, s1=1, sx=0.0267, sy=0.00034, sz=0.011)
@var Transforms.NTF: Transform(name='NTF', tx=-168, ty=-60, tz=320, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.OSGB36: Transform(name='OSGB36', tx=-446.448, ty=125.157, tz=-542.06, rx=-0, ry=-0, rz=-0, s=20.4894, s1=1.00002, sx=-0.1502, sy=-0.247, sz=-0.8421)
@var Transforms.TokyoJapan: Transform(name='TokyoJapan', tx=148, ty=-507, tz=-685, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.WGS72: Transform(name='WGS72', tx=0, ty=0, tz=-4.5, rx=0, ry=0, rz=0, s=-0.22, s1=1, sx=0, sy=0, sz=0.554)
@var Transforms.WGS84: Transform(name='WGS84', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
'''

# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .fmath.py, .utily.py
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.basics import EPS, EPS1, PI2, PI_2, R_M, \
                             property_doc_, property_RO, _xinstanceof
from pygeodesy.errors import _AssertionError, _ValueError
from pygeodesy.fmath import _2_3rd, cbrt, cbrt2, fdot, fpowers, Fsum, \
                             fsum_, hypot1, hypot2, sqrt3  # PYCHOK _2_3rd
from pygeodesy.interns import _COMMA_SPACE_, _ellipsoid_, _lat0_, _lat1_, \
                              _lon0_, _lon1_, _n_a_, _name_, NN, \
                              _transform_, _UNDERSCORE_, _x_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import Distance2Tuple, _NamedEnum, _NamedEnumItem, \
                           _NamedTuple, Vector3Tuple
from pygeodesy.streprs import _e, instr, _Fmt, fstr
from pygeodesy.units import Bearing_, Lam_, Lat, Phi, Phi_, \
                            Radius, Radius_, Scalar
from pygeodesy.utily import degrees360, degrees2m, m2degrees, \
                            m2km, m2NM, m2SM

from math import atan, atan2, atanh, copysign, cos, exp, hypot, \
                 radians, sin, sinh, sqrt

R_M  = Radius(R_M,       name='R_M')  #: Mean (spherical) earth radius (C{meter}).
R_MA = Radius(6378137.0, name='R_MA')  #: Major (equatorial) earth radius (C{meter}) WGS84, EPSG:3785.
R_MB = Radius(6356752.0, name='R_MB')  #: Minor (polar) earth radius (C{meter}) WGS84, EPSG:3785.
R_KM = Radius(m2km(R_M), name='R_KM')  #: Mean (spherical) earth radius (C{KM}, kilo meter).
R_NM = Radius(m2NM(R_M), name='R_NM')  #: Mean (spherical) earth radius (C{NM}, nautical miles).
R_SM = Radius(m2SM(R_M), name='R_SM')  #: Mean (spherical) earth radius (C{SM}, statute miles).
# See <https://www.EdWilliams.org/avform.htm>,
# <https://www.DTIC.mil/dtic/tr/fulltext/u2/a216843.pdf> and
# <https://GitHub.com/NASA/MultiDop/blob/master/src/share/man/man3/geog_lib.3>
# based on International Standard Nautical Mile of 1,852 meter (1' latitude)
R_FM = Radius(6371000.0,       name='R_FM')  #: Former FAI Sphere earth radius (C{meter}).
R_VM = Radius(6366707.0194937, name='R_VM')  #: Aviation/Navigation earth radius (C{meter}).
# R_ = Radius(6372797.560856,  name='R_')   #: XXX some other earth radius???

__all__ = _ALL_LAZY.datum
__version__ = '20.08.04'

_Flts = {}               # floats cache, deleted below
_TOL  = sqrt(EPS * 0.1)  # for Ellipsoid.estauf, imported by .ups


def _4Ecef(this, Ecef):
    '''Return an ECEF converter.
    '''
    from pygeodesy.ecef import EcefKarney, EcefVeness, EcefYou

    if Ecef is None:
        Ecef = EcefKarney
    else:
        _xinstanceof(EcefKarney, EcefVeness, EcefYou, Ecef=Ecef)
    return Ecef(this, name=this.name)  # datum or ellipsoid


def _flt(f):
    '''(INTERNAL) cache initial C{float}s.
    '''
    f = float(f)
    return _Flts.setdefault(f, f)  # PYCHOK del _Flts


class Curvature2Tuple(_NamedTuple):
    '''2-Tuple C{(meridional, prime_vertical)} of radii of curvature,
       both in C{meter}.
    '''
    _Names_ = ('meridional', 'prime_vertical')


class Ellipsoid(_NamedEnumItem):
    '''Ellipsoid with major and minor radius, semi-axis, (inverse)
       flattening and other pre-computed, frequently used attributes.
    '''
    _a  = 0  #: Major (equatorial) radius, semi-axis (C{meter}).
    _b  = 0  #: Minor (polar) radius, semi-axis (C{meter}): a * (f - 1) / f.
    _f  = 0  #: Flattening: (a - b) / a
    _f_ = 0  #: Inverse flattening: a / (a - b) = 1 / f

    _f2 = None  #: 2nd Flattening: (a - b) / b  # un-/rarely used
    _n  = None  #: 3rd Flattening: f / (2 - f) = (a - b) / (a + b)  # for .A and utm

    _a2  = None  #: a**2
    _a2_ = None  #: (1 / a**2)  # for ellipsiodalNvector.Cartesian.toNvector
    _a_b = None  #: (a / b) = 1 / (1 - f)  # for ellipsoidalNvector.Nvector.toCartesian
    _b2  = None  #: b**2
    _b_a = None  #: (b / a) = 1 - f  # for formy

    # curvatures <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>
    _a2_b = None  #: (INTERNAL) Meridional radius of curvature at poles: a**2 / b (C{meter})
    _b2_a = None  #: (INTERNAL) Meridional radius of curvature at equator: b**2 / a (C{meter})

    # eccentricities
    _e    = None  #: 1st Eccentricity: sqrt(1 - (b / a)**2))  # for utm
    _es_c = None  #: M{(1 - f) * exp(es_atanh(1))}
    _e2   = None  #: 1st Eccentricity squared: f * (2 - f) = 1 - (b / a)**2
    _e4   = None  #: e2**2 = e**4 for ellipsoidalNvector.Cartesian.toNvector, ecef
    _e12  = None  #: 1 - e2 = (1 - f)**2 for ellipsoidalNvector.Cartesian.toNvector, ecef utm
    _e22  = None  #: 2nd Eccentricity squared: e2 / (1 - e2) = (a / b)**2 - 1

    # fixed earth radii from <https://WikiPedia.org/wiki/Earth_radius>
    _R1 = None  #: (INTERNAL) Mean earth radius: (2 * a + b) / 3 per IUGG definition (C{meter})
    _R2 = None  #: (INTERNAL) Authalic radius: sqrt((a**2 + b**2 * atanh(e) / e) / 2) (C{meter})
#   _c  = None  #: (INTERNAL) Authalic radius: equ (60) in Karney's "Algorithms for Geodesics"
    _R3 = None  #: (INTERNAL) Volumetric radius: (a * a * b)**1/3 (C{meter})
    _Rr = None  #: (INTERNAL) Rectifying radius: ((a**3/2 + b**3/2) / 2)**2/3 (C{meter})
    _Rs = None  #: (INTERNAL) Mean earth radius: sqrt(a * b) (C{meter})

    _ab_90  = None  #: (a - b) / 90  # for .Rlat below
    _area   = None  #: (INTERNAL) Surface area: 4 * PI * R2**2
    _volume = None  #: (INTERNAL) Volume: 4 / 3 * PI * a**2 * b

    _A       = None  #: (INTERNAL) Meridional radius
    _AlphaKs = None  #: (INTERNAL) Up to 8th-order Krüger Alpha series
    _BetaKs  = None  #: (INTERNAL) Up to 8th-order Krüger Beta series
    _KsOrder = 8     #: (INTERNAL) Krüger series order (4, 6 or 8)
    _Mabcd   = None  #: (INTERNAL) OSGB meridional coefficients

    _geodesic = None  #: (INTERNAL) Cached C{karney._wrapped_.Geodesic} instance
    _Math     = None  #: (INTERNAL) Cached C{geographiclib.geomath.Math} module

    def __init__(self, a, b, f_, name=NN):
        '''New L{Ellipsoid}.

           @arg a: Semi-major, equatorial axis (C{meter}).
           @arg b: Semi-minor, polar axis (C{meter}).
           @arg f_: Inverse flattening: a / (a - b) (C{float} >>> 1.0).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}}, B{C{b}} or B{C{f_}}.
        '''
        self._a = a = _flt(a)  # major half-axis in meter
        if b:  # get f_ from a and b if not spherical
            self._b = b = _flt(b)
            if not f_ and a > b:
                f_ = a / (a - b)
        else:  # get b from a and f_, minor half-axis in meter
            self._b = b = _flt(a * (f_ - 1) / _flt(f_))

        if f_ > 0 and a > b > 0:
            self._f_ = f_ = _flt(f_)  # inverse flattening
            self._f  = f  = _flt(1 / f_)  # flattening
        elif a > 0:  # sphere
            self._b = b = self._a2b = self._b2a = a
            self._f2 = self._n = f_ = f = 0
            self._a_b = 1
            self._R1 = self._R2 = self._R3 = self._Rr = self._Rs = a
        else:
            raise _ValueError(instr(self, a, b, f_, name=name))

        d = a - b
        # some sanity checks to catch mistakes
        if d < 0 or min(a, b) < 1:
            raise _AssertionError('%s: %s=%0.9g vs %s=%0.9g' % (name,
                                  'a', a, 'b', b))
        t = d / a
        if abs(f - t) > 1e-8:
            raise _AssertionError('%s: %s=%.9g vs %s=%.9g' % (name,
                                  'f', f, '(a-b)/a', t))
        t = (1 - self.f)**2
        if abs(self.e12 - t) > EPS:
            raise _AssertionError('%s: %s=%.9g vs %s=%.9g' % (name,
                                  'e12', self.e12, '(1-f)**2', t))

        self._register(Ellipsoids, name)

    def __eq__(self, other):
        '''Compare this and an other ellipsoid.

           @arg other: The other ellipsoid (L{Ellipsoid}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Ellipsoid) and
                                 self.a == other.a and
                                (self.b == other.b or
                                 self.f == other.f))

    def _Kseries(self, *AB8Ks):
        '''(INTERNAL) Compute the 4-, 6- or 8-th order Krüger Alpha
           or Beta series coefficients per Karney 2011, 'Transverse
           Mercator with an accuracy of a few nanometers', U{page 7,
           equations 35 and 36<https://Arxiv.org/pdf/1002.1417v3.pdf>}.

           @arg AB8Ks: 8-Tuple of 8-th order Krüger Alpha or Beta series
                       coefficient tuples.

           @return: Krüger series coefficients (C{.KsOrder}-tuple).

           @see: The 30-th order U{TMseries30
                 <https://GeographicLib.SourceForge.io/html/tmseries30.html>}.
        '''
        k = self.KsOrder
        ns = fpowers(self.n, k)
        return tuple(fdot(AB8Ks[i][:k-i], *ns[i:]) for i in range(k))

    @property_RO
    def a(self):
        '''Get the major (equatorial) radius, semi-axis (C{meter}).
        '''
        return self._a

    @property_RO
    def a_b(self):
        '''Get ratio M{a / b} (C{float}).
        '''
        if self._a_b is None:
            self._a_b = self.a / self.b
        return self._a_b

    @property_RO
    def a2(self):
        '''Get the major radius I{squared} (C{float}), M{a**2}.
        '''
        if self._a2 is None:
            self._a2 = self.a**2
        return self._a2

    @property_RO
    def a2_(self):
        '''Get the inverse of the major radius I{squared} (C{float}), M{1 / a**2}.
        '''
        if self._a2_ is None:
            self._a2_ = 1 / self.a2
        return self._a2_  # (1 / a**2)

    @property_RO
    def a2_b(self):
        '''Get the polar meridional radius of curvature (C{meter}), M{a**2 / b}, see C{rocPolar}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
                 and U{Moritz, H. (1980), Geodetic Reference System 1980
                 <https://WikiPedia.org/wiki/Earth_radius#cite_note-Moritz-2>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of
                  curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        if self._a2_b is None:
            self._a2_b = self.a2 / self.b
        return self._a2_b

    @property_RO
    def area(self):
        '''Get the ellipsoid's surface area (C{meter**2}), M{4 * PI * R2**2}.
        '''
        if self._area is None:
            self._area = 2 * PI2 * self.R2**2
        return self._area

    @property_RO
    def A(self):
        '''Get the UTM meridional radius (C{meter}).
        '''
        if self._A is None:
            n = self.n
            # <https://GeographicLib.SourceForge.io/html/transversemercator.html>
            self._A = self.a / (1 + n) * (fsum_(65536, 16384 * n**2,
                                                        1024 * n**4,
                                                         256 * n**6,
                                                         100 * n**8,
                                                          49 * n**10) / 65536)
            # <https://www.MyGeodesy.id.AU/documents/Karney-Krueger%20equations.pdf>
            # self._A = self.a / (1 + n) * (fhorner(n**2, 16384, 4096, 256, 64, 25) / 16384)
        return self._A

    @property_RO
    def AlphaKs(self):
        '''Get the U{Krüger Alpha series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}-tuple).
        '''
        if self._AlphaKs is None:
            self._AlphaKs = self._Kseries(  # XXX int/int quotients may require  from __future__ import division
                # n    n**2   n**3      n**4         n**5            n**6                 n**7                     n**8
                (1/2, -2/3,   5/16,    41/180,    -127/288,       7891/37800,         72161/387072,        -18975107/50803200),
                     (13/48, -3/5,    557/1440,    281/630,   -1983433/1935360,       13769/28800,         148003883/174182400),      # PYCHOK unaligned
                            (61/240, -103/140,   15061/26880,   167603/181440,    -67102379/29030400,       79682431/79833600),       # PYCHOK unaligned
                                   (49561/161280, -179/168,    6601661/7257600,       97445/49896,      -40176129013/7664025600),     # PYCHOK unaligned
                                                (34729/80640, -3418889/1995840,    14644087/9123840,      2605413599/622702080),      # PYCHOK unaligned
                                                            (212378941/319334400, -30705481/10378368,   175214326799/58118860800),    # PYCHOK unaligned
                                                                                (1522256789/1383782400, -16759934899/3113510400),     # PYCHOK unaligned
                                                                                                      (1424729850961/743921418240,))  # PYCHOK unaligned
        return self._AlphaKs

    @property_RO
    def b(self):
        '''Get the minor (polar) radius, semi-axis (C{meter}).
        '''
        return self._b

    @property_RO
    def b_a(self):
        '''Get ratio M{b / a} (C{float}), M{1 - f}.
        '''
        if self._b_a is None:
            self._b_a = 1 - self.f
        return self._b_a

    @property_RO
    def b2(self):
        '''Get the minor radius I{squared} (C{float}), M{b**2}.
        '''
        if self._b2 is None:
            self._b2 = self.b**2
        return self._b2

    @property_RO
    def b2_a(self):
        '''Get the equatorial meridional radius of curvature (C{meter}), M{b**2 / a}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        if self._b2_a is None:
            self._b2_a = self.b2 / self.a
        return self._b2_a

    @property_RO
    def BetaKs(self):
        '''Get the U{Krüger Beta series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}-tuple).
        '''
        if self._BetaKs is None:
            self._BetaKs = self._Kseries(  # XXX int/int quotients may require  from __future__ import division
                # n    n**2  n**3     n**4        n**5            n**6                 n**7                   n**8
                (1/2, -2/3, 37/96,   -1/360,    -81/512,      96199/604800,     -5406467/38707200,      7944359/67737600),
                      (1/48, 1/15, -437/1440,    46/105,   -1118711/3870720,       51841/1209600,      24749483/348364800),       # PYCHOK unaligned
                           (17/480, -37/840,   -209/4480,      5569/90720,       9261899/58060800,     -6457463/17740800),        # PYCHOK unaligned
                                  (4397/161280, -11/504,    -830251/7257600,      466511/2494800,     324154477/7664025600),      # PYCHOK unaligned
                                              (4583/161280, -108847/3991680,    -8005831/63866880,     22894433/124540416),       # PYCHOK unaligned
                                                          (20648693/638668800, -16363163/518918400, -2204645983/12915302400),     # PYCHOK unaligned
                                                                              (219941297/5535129600, -497323811/12454041600),     # PYCHOK unaligned
                                                                                                  (191773887257/3719607091200,))  # PYCHOK unaligned
        return self._BetaKs

    @property_RO
    def c(self):
        '''Get the authalic earth radius (C{meter}), see C{R2}.

           @note: Symbol C{c} in U{equation 60
                  <https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of
                  curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>},
                  see C{rocPolar}.
        '''
        return self.R2 if self._R2 is None else self._R2

    def degrees2m(self, deg, lat=0):
        '''Convert angle to distance along the equator or along
           a parallel at an other latitude.

           @arg deg: Angle (C{degrees}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{deg}} or B{C{lat}}.
        '''
        return degrees2m(deg, radius=self.a, lat=lat)

    def distance2(self, lat0, lon0, lat1, lon1):
        '''Approximate the distance and (initial) bearing between two
           points based on the U{local, flat earth approximation
           <https://www.EdWilliams.org/avform.htm#flat>} aka U{Hubeny
           <https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           Suitable only for distances of several hundred Km or Miles
           and only between points not near-polar.

           @arg lat0: From latitude (C{degrees}).
           @arg lon0: From longitude (C{degrees}).
           @arg lat1: To latitude (C{degrees}).
           @arg lon1: To longitude (C{degrees}).

           @return: A L{Distance2Tuple}C{(distance, initial)} with C{distance}
                    in same units as this ellipsoid's axes.

           @note: The meridional and prime_vertical radii of curvature
                  are taken and scaled I{at the initial latitude}.

           @see: Function L{flatLocal}/L{hubeny}.
        '''
        phi0 = Phi_(lat0, name=_lat0_)
        m, n = self.roc2_(phi0, scaled=True)
        m *= Phi_(lat1, name=_lat1_) - phi0
        n *= Lam_(lon1, name=_lon1_) - Lam_(lon0, name=_lon0_)
        return Distance2Tuple(hypot(m, n), degrees360(atan2(n, m)))

    @property_RO
    def e(self):
        '''Get the (1st) Eccentricity (C{float}), M{sqrt(1 - (b / a)**2))}.
        '''
        if self._e is None:
            self._e = sqrt(abs(self.e2))
        return self._e

    @property_RO
    def e12(self):
        '''Get M{1 - e**2} (C{float}).
        '''
        if self._e12 is None:
            self._e12 = 1 - self.e2
        return self._e12  # 1 - e2

    @property_RO
    def e2(self):
        '''Get the (1st) Eccentricity I{squared} (C{float}), M{f * (2 - f) == 1 - (b / a)**2}.
        '''
        if self._e2 is None:
            self._e2 = e2 = self.f * (2 - self.f)
            t = 1 - (self.b / self.a)**2
            if abs(e2 - t) > 1e-9:
                raise _AssertionError('%s: %s=%.9g vs %s=%.9g' % (self.name,
                                      'e2', e2, '1-(b/a)**2', t))
        return self._e2

    @property_RO
    def e22(self):
        '''Get the 2nd Eccentricity I{squared} (C{float}), M{e2 / (1 - e2) == (a / b)**2 - 1}.
        '''
        if self._e22 is None:
            self._e22 = e22 = self.e2 / (1 - self.e2)
            t = self.a_b**2 - 1
            if abs(e22 - t) > 1e-9:
                raise _AssertionError('%s: %s=%.9g vs %s=%.9g' % (self.name,
                                      'e22', e22, '(a/b)**2-1', t))
        return self._e22

    @property_RO
    def e4(self):
        '''Get the (1st) Eccentricity to 4th power (C{float}), M{e**4 == e2**2}.
        '''
        if self._e4 is None:
            self._e4 = self.e2**2
        return self._e4

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{ellipsoid} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    def e2s(self, s):
        '''Compute norm M{sqrt(1 - e2 * s**2)}.

           @arg s: S value (C{scalar}).

           @return: Norm (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        return sqrt(self.e2s2(s))

    def e2s2(self, s):
        '''Compute M{1 - e2 * s**2}.

           @arg s: S value (C{scalar}).

           @return: Result (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        try:
            r = 1 - self.e2 * Scalar(s, name='s')**2
            if r < 0:
                raise ValueError
        except (TypeError, ValueError) as x:
            raise _ValueError(self._dot_('e2s2'), s, txt=str(x))
        return r

    @property_RO
    def es(self):
        '''Get the (1st) Eccentricity I{signed} (C{float}).
        '''
        # note, self.e is always non-negative
        return copysign(self.e, self.f)  # see .ups.py

    def es_atanh(self, x):
        '''Compute M{es * atanh(es * x)} where I{es} is the I{signed}
           1st Eccentricity.

           @raise ValueError: Invalid B{C{x}}.

           @see: Function U{Math::eatanhe<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        # note, self.e is always non-negative
        if self.f > 0:
            r = self.e * atanh(self.e * Scalar(x, name=_x_))
        elif self.f < 0:  # PYCHOK no cover
            r = self.e * atan(-self.e * Scalar(x, name=_x_))
        else:
            r = 0
        return r

    @property_RO
    def es_c(self):
        '''Get M{(1 - f) * exp(es_atanh(1))} (C{float}), M{b_a * exp(es_atanh(1))}.
        '''
        if self._es_c is None:
            self._es_c = self.b_a * exp(self.es_atanh(1.0))
        return self._es_c

    def es_tauf(self, taup):
        '''Compute U{Karney's<https://ArXiv.org/abs/1002.1417>}
           equations (19), (20) and (21).

           @see: Function U{Math::tauf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        T_  = Scalar(taup, name='taup')
        tol = max(abs(T_), 1) * _TOL
        e2m = 1 - abs(self.e2)  # == self.e**2
        t = T_ / e2m
        T = Fsum(t)
        for _ in range(5):
            a = self.es_taupf(t)
            d = (T_ - a) * (1 + e2m * t**2) / (e2m * hypot1(t) * hypot1(a))
            t, d = T.fsum2_(d)
            if abs(d) < tol:
                break
        return t

    def es_taupf(self, tau):
        '''Compute U{Karney's<https://ArXiv.org/abs/1002.1417>}
           equations (7), (8) and (9).

           @see: Function U{Math::taupf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        T = Scalar(tau, name='tau')
        t = hypot1(T)
        s = sinh(self.es_atanh(T / t))
        return hypot1(s) * T - s * t

    @property_RO
    def f(self):
        '''Get the Flattening (C{float}), M{(a - b) / a}.
        '''
        return self._f

    @property_RO
    def f_(self):
        '''Get the Inverse flattening (C{float}), M{1 / f} or M{a / (a - b)}.
        '''
        return self._f_

    @property_RO
    def f2(self):
        '''Get the 2nd Flattening (C{float}), M{(a - b) / b}.
        '''
        if self._f2 is None:
            self._f2 = (self.a - self.b) / self.b
        return self._f2

    @property_RO
    def geodesic(self):
        '''Get this ellipsoid's I{wrapped} U{Karney Geodesic
           <https://GeographicLib.SourceForge.io/html/python/code.html>},
           provided the U{GeographicLib
           <https://PyPI.org/project/geographiclib>} package is installed.
        '''
        if self._geodesic is None:
            # if not self.isEllipsoidal:
            #     raise _IsnotError('ellipsoidal', ellipsoid=self)
            from pygeodesy.karney import _wrapped
            self._geodesic = _wrapped.Geodesic(self.a, self.f)
        return self._geodesic

    @property_RO
    def _geodesic_Math2(self):
        '''(INTERNAL) Get this ellipsoid's C{Geodesic} and I{Karney}'s
           C{Math} class.
        '''
        if Ellipsoid._Math is None:
            from pygeodesy.karney import _wrapped
            Ellipsoid._Math = _wrapped.Math
        return self.geodesic, Ellipsoid._Math

    def _hubeny2_(self, phi2, phi1, lam21):
        '''(INTERNAL) like function L{flatLocal_}/L{hubeny_} but
           returning the I{angular} distance in C{radians squared}.
        '''
        m, n = self.roc2_((phi2 + phi1) * 0.5, scaled=True)
        return hypot2(m * (phi2 - phi1), n * lam21) * self.a2_

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this model is ellipsoidal (C{bool}).
        '''
        return self.a > self.R1 > self.b

    @property_RO
    def isSpherical(self):
        '''Check whether this model is spherical (C{bool}).
        '''
        return self.a == self.R1 == self.b

    @property_doc_(''' the Krüger series' order (C{int}).''')
    def KsOrder(self):
        '''Get the Krüger series order (C{int} 4, 6 or 8).
        '''
        return self._KsOrder

    @KsOrder.setter  # PYCHOK setter!
    def KsOrder(self, order):
        '''Set the Krüger series' order.

           @arg order: New Krüger series' order (C{int} 4, 6 or 8).

           @raise ValueError: Invalid B{C{order}}.
        '''
        if order not in (4, 6, 8):
            raise _ValueError(order=order)
        if order != self._KsOrder:
            if self._AlphaKs:
                self._AlphaKs = None
            if self._BetaKs:
                self._BetaKs = None
            self._KsOrder = order

    def m2degrees(self, meter, lat=0):
        '''Convert distance to angle along equator or along
           a parallel at an other latitude.

           @arg meter: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{meter}} or B{C{lat}}.
       '''
        return m2degrees(meter, radius=self.a, lat=lat)

    @property_RO
    def Mabcd(self):
        '''Get the OSGR meridional coefficients (C{4-Tuple}), C{Airy130} only.
        '''
        if self._Mabcd is None:
            n, n2, n3 = fpowers(self.n, 3)  # PYCHOK false!
            self._Mabcd = (fdot((1, n, n2, n3), 4, 4, 5, 5) / 4,
                           fdot(   (n, n2, n3), 24, 24, 21) / 8,
                           fdot(      (n2, n3), 15, 15) / 8,
                                      35 * n3 / 24)
        return self._Mabcd

    @property_RO
    def n(self):
        '''Get the 3rd Flattening (C{float}), M{f / (2 - f) == (a - b) / (a + b)}.
        '''
        if self._n is None:
            self._n = n = self.f / (2 - self.f)
            t = (self.a - self.b) / (self.a + self.b)
            if abs(n - t) > 1e-9:
                raise _AssertionError('%s: %s=%.9g vs %s=%.9g' % (self.name,
                                      'n', n, '(a-b)/(a+b)', t))
        return self._n

    @property_RO
    def R1(self):
        '''Get the mean earth radius per IUGG (C{meter}), M{(2 * a + b) / 3}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        if self._R1 is None:
            self._R1 = Radius((self.a * 2 + self.b) / 3, name=Ellipsoid.R1.name)
        return self._R1

    @property_RO
    def R2(self):
        '''Get the authalic earth radius (C{meter}), M{sqrt((a**2 + b**2 * atanh(e) / e) / 2)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>} and
                 U{c<https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>}.
        '''
        if self._R2 is None:
            if self.e2 > 0:
                r = atanh(self.e) / self.e
            elif self.e2 < 0:  # GeographicLib Ellipsoid Area
                r = atan(self.e) / self.e
            else:
                r = 1
            self._R2 = Radius(sqrt((self.a2 + self.b2 * r) * 0.5), name=Ellipsoid.R2.name)
        return self._R2

    @property_RO
    def R3(self):
        '''Get the volumetric earth radius (C{meter}), M{(a * a * b)**(1/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        if self._R3 is None:
            self._R3 = Radius(cbrt(self.a2 * self.b), name=Ellipsoid.R3.name)
        return self._R3

    def Rgeocentric(self, lat):
        '''Compute the geocentric earth radius at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Geocentric earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Geocentric Radius
                 <https://WikiPedia.org/wiki/Earth_radius#Geocentric_radius>}
        '''
        a2 = self.a2
        b2 = self.b2
        c2 = cos(Phi_(lat))**2
        s2 = 1 - c2
        return Radius(sqrt((a2**2 * c2 + b2**2 * s2) / (a2 * c2 + b2 * s2)),
                      name=Ellipsoid.Rgeocentric.__name__)

    @property_RO
    def Rr(self):
        '''Get the rectifying earth radius (C{meter}), M{((a**(3/2) + b**(3/2)) / 2)**(2/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        if self._Rr is None:
            self._Rr = Radius(cbrt2((sqrt3(self.a) + sqrt3(self.b)) * 0.5),
                              name=Ellipsoid.Rr.name)
        return self._Rr

    @property_RO
    def Rs(self):
        '''Get another mean earth radius (C{meter}), M{sqrt(a * b)}.
        '''
        if self._Rs is None:
            self._Rs = Radius(sqrt(self.a * self.b), name=Ellipsoid.Rs.name)
        return self._Rs

    def Rlat(self, lat):
        '''Approximate the earth radius at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Approximate earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.
        '''
        if self._ab_90 is None:
            self._ab_90 = (self.a - self.b) / 90.0
        # r = major - (major - minor) * |lat| / 90
        r = self.a
        if lat:
            r -= self._ab_90 * min(abs(Lat(lat, clip=0)), 90)
        return Radius(r, name=Ellipsoid.Rlat.__name__)

    def roc2(self, lat, scaled=False):
        '''Compute the meridional and prime-vertical, I{normal} radii of curvature
           at the given latitude.

           @arg lat: Latitude (C{degrees90}).
           @kwarg scaled: Scale prime_vertical by B{C{cos(phi)}} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Method C{.roc2_}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat), scaled=scaled)

    def roc2_(self, phi, scaled=False):
        '''Compute the meridional and prime-vertical, I{normal} radii of curvature
           at the given latitude.

           @arg phi: Latitude (C{radians}).
           @kwarg scaled: Scale prime_vertical by B{C{cos(phi)}} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{phi}}.

           @see: Method C{.roc2}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and meridional
                 and prime vertical U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        a = abs(Phi(phi))
        r = self.e2s2(sin(a) if a < PI_2 else 1)
        if r < EPS:
            m = n = 0  # PYCHOK e2s2 attr?
        elif r < EPS1:
            n = self.a / sqrt(r)
            m = n * self.e12 / r  # PYCHOK e2s2 attr?
        else:
            n = self.a
            m = n * self.e12
        if scaled:
            n *= cos(a) if a < PI_2 else 0
        return Curvature2Tuple(Radius(m, name=Ellipsoid.rocMeridional.__name__),
                               Radius(n, name=Ellipsoid.rocPrimeVertical.__name__))

    def rocBearing(self, lat, bearing):
        '''Compute the directional radius of curvature at the
           given latitude and compass direction.

           @arg lat: Latitude (C{degrees90}).
           @arg bearing: Direction (compass C{degrees360}).

           @return: Directional radius of curvature (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{bearing}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        c2 = cos(Bearing_(bearing))**2
        s2 = 1 - c2
        m, n = self.roc2_(Phi_(lat))
        if n < m:  # == n / (c2 * n / m + s2)
            c2 *= n / m
        elif m < n:  # == m / (c2 + s2 * m / n)
            s2 *= m / n
            n = m
        return Radius(n / (c2 + s2), name=Ellipsoid.rocBearing.__name__)  # == 1 / (c2 / m + s2 / n)

    def rocGauss(self, lat):
        '''Compute the Gaussian radius of curvature at the given latitude.

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
        c2 = cos(Phi_(lat))**2
        return Radius(self.b / (c2 + (1 - c2) * self.b2 / self.a2), name=Ellipsoid.rocGauss.__name__)

    def rocMean(self, lat):
        '''Compute the mean radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Mean radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        m, n = self.roc2_(Phi_(lat))
        return Radius(2 * m * n / (m + n), name=Ellipsoid.rocMean.__name__)  # == 2 / (1 / m + 1 / n)

    def rocMeridional(self, lat):
        '''Compute the meridional radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Meridional radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).meridional

    @property_RO
    def rocPolar(self):
        '''Get the polar radius of curvature (C{meter}), see C{a2_b}.
        '''
        return self.a2_b if self._a2_b is None else self._a2_b

    def rocPrimeVertical(self, lat):
        '''Compute the prime-vertical, I{normal} radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Prime-vertical radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).prime_vertical

    def toStr(self, prec=9):  # PYCHOK expected
        '''Return this ellipsoid as a text string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Ellipsoid attributes (C{str}).
        '''
        return self._instr(prec, 'a', 'b', 'f_', 'f', 'e', 'e2', 'e12', 'e22',
                                 'n', 'R1', 'R2', 'R3', 'Rr', 'Rs')

    @property_RO
    def volume(self):
        '''Get the ellipsoid's volume (C{meter**3}), M{4 / 3 * PI * a**2 * b}.
        '''
        if self._volume is None:
            self._volume = Scalar(PI2 * _2_3rd * self.a2 * self.b, name=Ellipsoid.volume.name)
        return self._volume


Ellipsoids = _NamedEnum('Ellipsoids', Ellipsoid)  #: Registered ellipsoids.
# <https://www.GNU.org/software/gama/manual/html_node/Supported-ellipsoids.html>
# <https://w3.Energistics.org/archive/Epicentre/Epicentre_v3.0/DataModel/
#         LogicalDictionary/StandardValues/ellipsoid.html>
# <https://kb.OSU.edu/dspace/handle/1811/77986>
Ellipsoids._assert(  # <https://WikiPedia.org/wiki/Earth_ellipsoid>
    Airy1830       = Ellipsoid(6377563.396, 6356256.909,       299.3249646,   'Airy1830'),
    AiryModified   = Ellipsoid(6377340.189, 6356034.448,       299.3249646,   'AiryModified'),
#   ANS            = Ellipsoid(6378160.0,   6356774.719,       298.25,        'ANS'),
    Australia1966  = Ellipsoid(6378160.0,   6356774.719,       298.25,        'Australia1966'),
#   Bessel1841     = Ellipsoid(6377397.155, 6356078.963,       299.152815351, 'Bessel1841'),
    Bessel1841     = Ellipsoid(6377397.155, 6356078.962818,    299.1528128,   'Bessel1841'),
    Clarke1866     = Ellipsoid(6378206.4,   6356583.8,         294.978698214, 'Clarke1866'),
    Clarke1880     = Ellipsoid(6378249.145, 6356514.86954978,  293.465,       'Clarke1880'),
    Clarke1880IGN  = Ellipsoid(6378249.2,   6356515.0,         293.466021294, 'Clarke1880IGN'),
    Clarke1880Mod  = Ellipsoid(6378249.145, 6356514.96582849,  293.4663,      'Clarke1880Mod'),
    CPM1799        = Ellipsoid(6375738.7,   6356671.92557493,  334.39,        'CPM1799'),  # Comm. des Poids et Mesures
    Delambre1810   = Ellipsoid(6376428.0,   6355957.92616372,  311.5,         'Delambre1810'),  # Belgium
    Engelis1985    = Ellipsoid(6378136.05,  6356751.32272154,  298.2566,      'Engelis1985'),
    Everest1969    = Ellipsoid(6377295.664, 6356094.667915,    300.8017,      'Everest1969'),
    Fisher1968     = Ellipsoid(6378150.0,   6356768.33724438,  298.3,         'Fisher1968'),
    GEM10C         = Ellipsoid(6378137.0,   6356752.31424783,  298.2572236,   'GEM10C'),
    GRS67          = Ellipsoid(6378160.0,   6356774.516,       298.247167427, 'GRS67'),  # Lucerne
    GRS80          = Ellipsoid(6378137.0,   6356752.314140347, 298.257222101, 'GRS80'),  # ITRS, ETRS89
    Helmert1906    = Ellipsoid(6378200.0,   6356818.16962789,  298.3,         'Helmert1906'),
    IERS1989       = Ellipsoid(6378136.0,   6356751.302,       298.257,       'IERS1989'),
    IERS1992TOPEX  = Ellipsoid(6378136.3,   6356751.61659215,  298.257223563, 'IERS1992TOPEX'),  # IERS/TOPEX/Poseidon/McCarthy
    IERS2003       = Ellipsoid(6378136.6,   6356751.85797165,  298.25642,     'IERS2003'),
    Intl1924       = Ellipsoid(6378388.0,   6356911.946,       297.0,         'Intl1924'),  # aka Hayford
    Intl1967       = Ellipsoid(6378157.5,   6356772.2,         298.24961539,  'Intl1967'),  # New Int'l
    Krassovski1940 = Ellipsoid(6378245.0,   6356863.01877305,  298.3,         'Krassovski1940'),  # spelling
    Krassowsky1940 = Ellipsoid(6378245.0,   6356863.01877305,  298.3,         'Krassowsky1940'),  # spelling
    Maupertuis1738 = Ellipsoid(6397300.0,   6363806.28272251,  191.0,         'Maupertuis1738'),  # France
    Mercury1960    = Ellipsoid(6378166.0,   6356784.28360711,  298.3,         'Mercury1960'),
    Mercury1968Mod = Ellipsoid(6378150.0,   6356768.33724438,  298.3,         'Mercury1968Mod'),
    NWL1965        = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'NWL1965'),  # Naval Weapons Lab.
    OSU86F         = Ellipsoid(6378136.2,   6356751.51693008,  298.2572236,   'OSU86F'),
    OSU91A         = Ellipsoid(6378136.3,   6356751.6165948,   298.2572236,   'OSU91A'),
#   Plessis1817    = Ellipsoid(6397523.0,   6355863.0,         153.56512242,  'Plessis1817'),  # XXX incorrect?
    Plessis1817    = Ellipsoid(6376523.0,   6355862.93325557,  308.64,        'Plessis1817'),  # XXX IGN France 1972
    SGS85          = Ellipsoid(6378136.0,   6356751.30156878,  298.257,       'SGS85'),  # Soviet Geodetic System
    SoAmerican1969 = Ellipsoid(6378160.0,   6356774.71919531,  298.25,        'SoAmerican1969'),  # South American
    Struve1860     = Ellipsoid(6378298.3,   6356657.14266956,  294.73,        'Struve1860'),
    WGS60          = Ellipsoid(6378165.0,   6356783.28695944,  298.3,         'WGS60'),
    WGS66          = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'WGS66'),
    WGS72          = Ellipsoid(6378135.0,   6356750.52,        298.26,        'WGS72'),
    WGS84          = Ellipsoid(6378137.0,   6356752.31425,     298.257223563, 'WGS84'),  # GPS
    Sphere         = Ellipsoid(R_M,         R_M,                 0.0,         'Sphere'),  # pseudo
    SphereAuthalic = Ellipsoid(R_FM,        R_FM,                0.0,         'SphereAuthalic'),  # pseudo
    SpherePopular  = Ellipsoid(R_MA,        R_MA,                0.0,         'SpherePopular'),  # EPSG:3857 Spheroid
)


def _r_s2(s):
    '''(INTERNAL) rotation in C{radians} and C{degree seconds}.
    '''
    return _flt(radians(s / 3600.0)), _flt(s)


class Transform(_NamedEnumItem):
    '''Helmert transformation.
    '''
    tx = 0  #: X translation (C{meter}).
    ty = 0  #: Y translation (C{meter}).
    tz = 0  #: Z translation (C{meter}).

    rx = 0  #: X rotation (C{radians}).
    ry = 0  #: Y rotation (C{radians}).
    rz = 0  #: Z rotation (C{radians}).

    s  = 0  #: Scale ppm (C{float}).
    s1 = 1  #: Scale + 1 (C{float}).

    sx = 0  #: X rotation (degree seconds).
    sy = 0  #: Y rotation (degree seconds).
    sz = 0  #: Z rotation (degree seconds).

    def __init__(self, name=NN, tx=0, ty=0, tz=0,
                                sx=0, sy=0, sz=0, s=0):
        '''New L{Transform}.

           @kwarg name: Optional, unique name (C{str}).
           @kwarg tx: Optional X translation (C{meter}).
           @kwarg ty: Optional Y translation (C{meter}).
           @kwarg tz: Optional Z translation (C{meter}).
           @kwarg s: Optional scale ppm (C{float}).
           @kwarg sx: Optional X rotation (C{degree seconds}).
           @kwarg sy: Optional Y rotation (C{degree seconds}).
           @kwarg sz: Optional Z rotation (C{degree seconds}).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        if tx:
            self.tx = _flt(tx)
        if ty:
            self.ty = _flt(ty)
        if tz:
            self.tz = _flt(tz)
        if sx:  # secs to rads
            self.rx, self.sx = _r_s2(sx)
        if sy:
            self.ry, self.sy = _r_s2(sy)
        if sz:
            self.rz, self.sz = _r_s2(sz)
        if s:
            self.s  = _flt(s)
            self.s1 = _flt(s * 1e-6 + 1)  # normalize ppm to (s + 1)

        self._register(Transforms, name)

    def __eq__(self, other):
        '''Compare this and an other transform.

           @arg other: The other transform (L{Transform}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Transform) and
                                 self.tx == other.tx and
                                 self.ty == other.ty and
                                 self.tz == other.tz and
                                 self.rx == other.rx and
                                 self.ry == other.ry and
                                 self.rz == other.rz and
                                 self.s  == other.s)

    def inverse(self, name=NN):
        '''Return the inverse of this transform.

           @kwarg name: Optional, unique name (C{str}).

           @return: Inverse (Transform).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        return Transform(name=name or (self.name + 'Inverse'),
                         tx=-self.tx, ty=-self.ty, tz=-self.tz,
                         sx=-self.sx, sy=-self.sy, sz=-self.sz, s=-self.s)

    def toStr(self, prec=5):  # PYCHOK expected
        '''Return this transform as a string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Transform attributes (C{str}).
        '''
        return self._instr(prec, 'tx', 'ty', 'tz',
                                 'rx', 'ry', 'rz', 's', 's1',
                                 'sx', 'sy', 'sz')

    def transform(self, x, y, z, inverse=False):
        '''Transform a (geocentric) Cartesian point, forward or inverse.

           @arg x: X coordinate (C{meter}).
           @arg y: Y coordinate (C{meter}).
           @arg z: Z coordinate (C{meter}).
           @kwarg inverse: Optional direction, forward or inverse (C{bool}).

           @return: A L{Vector3Tuple}C{(x, y, z)}, transformed.
        '''
        if inverse:
            _xyz = -1, -x, -y, -z
            _s1 = self.s1 - 2  # == -(1 - s * 1e-6)) == -(1 - (s1 - 1))
        else:
            _xyz =  1,  x,  y,  z
            _s1  = self.s1
        # x', y', z' = (.tx + x * .s1 - y * .rz + z * .ry,
        #               .ty + x * .rz + y * .s1 - z * .rx,
        #               .tz - x * .ry + y * .rx + z * .s1)
        r = Vector3Tuple(fdot(_xyz, self.tx,      _s1, -self.rz,  self.ry),
                         fdot(_xyz, self.ty,  self.rz,      _s1, -self.rx),
                         fdot(_xyz, self.tz, -self.ry,  self.rx,      _s1))
        return self._xnamed(r)


Transforms = _NamedEnum('Transforms', Transform)  #: Registered transforms.
# <https://WikiPedia.org/wiki/Helmert_transformation> from WGS84
Transforms._assert(
    BD72           = Transform('BD72', tx=106.868628, ty=-52.297783, tz=103.723893,
                     # <https://www.NGI.BE/FR/FR4-4.shtm> ETRS89 == WG84
                     # <https://GeoRepository.com/transformation_15929/BD72-to-WGS-84-3.html>
                                       sx=-0.33657,   sy= -0.456955, sz= -1.84218,
                                        s= 1.2727),
    Bessel1841     = Transform('Bessel1841', tx=-582.0,  ty=-105.0, tz=-414.0,
                                             sx=  -1.04, sy= -0.35, sz=   3.08,
                                              s=  -8.3),
    Clarke1866     = Transform('Clarke1866', tx=8, ty=-160, tz=-176),
    DHDN           = Transform('DHDN', tx=-591.28,  ty=-81.35,   tz=-396.39,
                                       sx=   1.477, sy= -0.0736, sz=  -1.458,
                                        s=  -9.82),  # Germany
    ED50           = Transform('ED50', tx=89.5, ty=93.8, tz=123.1,
                     # <https://GeoNet.ESRI.com/thread/36583> sz=-0.156
                     # <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal.js>
                     # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
                                                         sz=  0.156, s=-1.2),
    Irl1965        = Transform('Irl1965', tx=-482.530, ty=130.596, tz=-564.557,
                                          sx=   1.042, sy=  0.214, sz=   0.631,
                                           s=  -8.15),
    Irl1975        = Transform('Irl1975', tx=-482.530, ty=130.596, tz=-564.557,
                     # XXX rotation signs may be opposite, to be checked
                                          sx=  -1.042, sy= -0.214, sz=  -0.631,
                                           s=  -1.1),
    Krassovski1940 = Transform('Krassovski1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),  # spelling
    Krassowsky1940 = Transform('Krassowsky1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),  # spelling
    MGI            = Transform('MGI', tx=-577.326, ty=-90.129, tz=-463.920,
                                      sx=   5.137, sy=  1.474, sz=   5.297,
                                       s=  -2.423),  # Austria
    NAD27          = Transform('NAD27', tx=8, ty=-160, tz=-176),
    NAD83          = Transform('NAD83', tx= 1.004,  ty=-1.910,   tz=-0.515,
                                        sx= 0.0267, sy= 0.00034, sz= 0.011,
                                         s=-0.0015),
    NTF            = Transform('NTF', tx=-168, ty= -60, tz=320),  # XXX verify
    OSGB36         = Transform('OSGB36', tx=-446.448,  ty=125.157,  tz=-542.060,
                                         sx=  -0.1502, sy= -0.2470, sz=  -0.8421,
                                          s=  20.4894),
    TokyoJapan     = Transform('TokyoJapan', tx=148, ty=-507, tz=-685),
    WGS72          = Transform('WGS72', tz=-4.5, sz=0.554, s=-0.22),
    WGS84          = Transform('WGS84'),  # unity
)


class Datum(_NamedEnumItem):
    '''Ellipsoid and transform parameters for an earth model.
    '''
    _ellipsoid = Ellipsoids.WGS84  #: (INTERNAL) Default ellipsoid (L{Ellipsoid}).
    _exactTM   = None              #: (INTERNAL) L{ExactTransverseMercator} projection.
    _transform = Transforms.WGS84  #: (INTERNAL) Default transform (L{Transform}).

    def __init__(self, ellipsoid, transform=None, name=NN):
        '''New L{Datum}.

           @arg ellipsoid: The ellipsoid (L{Ellipsoid}).
           @kwarg transform: Optional transform (L{Transform}).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Datum with that B{C{name}} already exists.

           @raise TypeError: If B{C{ellipsoid}} is not an L{Ellipsoid}
                             or B{C{transform}} is not a L{Transform}.
        '''
        self._ellipsoid = ellipsoid or Datum._ellipsoid
        _xinstanceof(Ellipsoid, ellipsoid=self.ellipsoid)

        self._transform = transform or Datum._transform
        _xinstanceof(Transform, transform=self.transform)

        self._register(Datums, name or self.transform.name or self.ellipsoid.name)

    def __eq__(self, other):
        '''Compare this and an other datum.

           @arg other: The other datum (L{Datum}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Datum) and
                                 self.ellipsoid == other.ellipsoid and
                                 self.transform == other.transform)

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{datum} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    @property_RO
    def ellipsoid(self):
        '''Get this datum's ellipsoid (L{Ellipsoid}).
        '''
        return self._ellipsoid

    @property_RO
    def exactTM(self):
        '''Get the C{ExactTM} projection (L{ExactTransverseMercator}).
        '''
        if self._exactTM is None:
            from pygeodesy.etm import ExactTransverseMercator
            self._exactTM = ExactTransverseMercator(datum=self)
        return self._exactTM

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this datum is ellipsoidal (C{bool}).
        '''
        return self._ellipsoid.isEllipsoidal

    @property_RO
    def isSpherical(self):
        '''Check whether this datum is spherical (C{bool}).
        '''
        return self._ellipsoid.isSpherical

    def toStr(self, **unused):  # PYCHOK expected
        '''Return this datum as a string.

           @return: Datum attributes (C{str}).
        '''
        t = ['%s=%r' % (_name_, self.named)]
        for a in (_ellipsoid_, _transform_):
            v = getattr(self, a)
            t.append('%s=%ss.%s' % (a, v.classname, v.name))
        return _COMMA_SPACE_.join(t)

    @property_RO
    def transform(self):
        '''Get this datum's transform (L{Transform}).
        '''
        return self._transform


def _spherical_datum(radius, name=NN):
    '''(INTERNAL) Create a spherical ellipsoid and datum.
    '''
    n = _UNDERSCORE_ + name
    r = Radius_(radius)
    E = Ellipsoid(r, r, 0, name=n)
    return Datum(E, Transforms.WGS84, name=n)


Datums = _NamedEnum('Datums', Datum)      #: Registered datums.
# Datums with associated ellipsoid and Helmert transform parameters
# to convert from WGS84 into the given datum.  More are available at
# <https://Earth-Info.NGA.mil/GandG/coordsys/datums/NATO_DT.pdf> and
# <XXX://www.FieldenMaps.info/cconv/web/cconv_params.js>.
Datums._assert(
    # Belgian Datum 1972, based on Hayford ellipsoid.
    # <https://NL.WikiPedia.org/wiki/Belgian_Datum_1972>
    # <https://SpatialReference.org/ref/sr-org/7718/html/>
    BD72           = Datum(Ellipsoids.Intl1924, Transforms.BD72),

    # Germany <https://WikiPedia.org/wiki/Bessel-Ellipsoid>
    #         <https://WikiPedia.org/wiki/Helmert_transformation>
    DHDN           = Datum(Ellipsoids.Bessel1841, Transforms.DHDN),

    # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
    ED50           = Datum(Ellipsoids.Intl1924, Transforms.ED50),

    # Australia <https://ICSM.Gov.AU/datum/gda2020-and-gda94-technical-manuals>
#   ADG66          = Datum(Ellipsoids.ANS, Transforms.WGS84, name='ADG66'),  # XXX Transform?
#   ADG84          = Datum(Ellipsoids.ANS, Transforms.WGS84, name='ADG84'),  # XXX Transform?
#   GDA94          = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GDA94'),
    GDA2020        = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GDA2020'),  # XXX Transform?

    # <https://WikiPedia.org/wiki/GRS_80>
    GRS80          = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GRS80'),

    # <https://OSI.IE/OSI/media/OSI/Content/Publications/transformations_booklet.pdf>
    Irl1975        = Datum(Ellipsoids.AiryModified, Transforms.Irl1975),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Krassovski1940 = Datum(Ellipsoids.Krassovski1940, Transforms.Krassovski1940),  # XXX spelling?
    Krassowsky1940 = Datum(Ellipsoids.Krassowsky1940, Transforms.Krassowsky1940),  # XXX spelling?

    # Austria <https://DE.WikiPedia.org/wiki/Datum_Austria>
    MGI            = Datum(Ellipsoids.Bessel1841, Transforms.MGI),

    # <https://WikiPedia.org/wiki/Helmert_transformation>
    NAD27          = Datum(Ellipsoids.Clarke1866, Transforms.NAD27),

    # NAD83 (2009) == WGS84 - <https://www.UVM.edu/giv/resources/WGS84_NAD83.pdf>
    # (If you *really* must convert WGS84<->NAD83, you need more than this!)
    NAD83          = Datum(Ellipsoids.GRS80, Transforms.NAD83),

    #  Nouvelle Triangulation Francaise (Paris)  XXX verify
    NTF            = Datum(Ellipsoids.Clarke1880IGN, Transforms.NTF),

    # <https://www.OrdnanceSurvey.co.UK/docs/support/guide-coordinate-systems-great-britain.pdf>
    OSGB36         = Datum(Ellipsoids.Airy1830, Transforms.OSGB36),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Potsdam        = Datum(Ellipsoids.Bessel1841, Transforms.Bessel1841, name='Potsdam'),

    # XXX psuedo-ellipsoids for spherical LatLon
    Sphere         = Datum(Ellipsoids.Sphere, Transforms.WGS84, name='Sphere'),

    # <https://www.GeoCachingToolbox.com?page=datumEllipsoidDetails>
    TokyoJapan     = Datum(Ellipsoids.Bessel1841, Transforms.TokyoJapan),

    # <https://www.ICAO.int/safety/pbn/documentation/eurocontrol/eurocontrol%20wgs%2084%20implementation%20manual.pdf>
    WGS72          = Datum(Ellipsoids.WGS72, Transforms.WGS72),

    WGS84          = Datum(Ellipsoids.WGS84, Transforms.WGS84),
)

__all__ += _ALL_DOCS(Curvature2Tuple)

if __name__ == '__main__':

    for E in (Datums.WGS84.ellipsoid, Datums.NAD83.ellipsoid,
              Ellipsoids.Sphere, Ellipsoids.SpherePopular):
        if E.f:
            f_ = 'f_=1/%.10F' % (1 / E.f,)
        else:
            f_ = 'f_=' + _n_a_
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = (E.toStr(prec=10),
            'A=%r, e=%s, %s, n=%s(%s)' % (E.A, fstr(E.e, prec=13, fmt=_e),
                                          f_,  fstr(E.n, prec=13, fmt=_Fmt),
                                               fstr(e,   prec=3,  fmt=_e),),
            '%s=(%s)'   % (Ellipsoid.AlphaKs.name, fstr(E.AlphaKs, prec=20),),
            '%s= (%s)'  % (Ellipsoid.BetaKs.name,  fstr(E.BetaKs,  prec=20),),
            '%s= %s'    % ('KsOrder',                   E.KsOrder),
            '%s=  (%s)' % (Ellipsoid.Mabcd.name,   fstr(E.Mabcd,   prec=20),))
        print('\nEllipsoid.%s: %s' % (E.name, ',\n    '.join(t)))

    # __doc__ of this file
    for e in (Datums, Ellipsoids, Transforms):
        t = [NN] + repr(e).split('\n')
        print('\n@var '.join(i.strip(',') for i in t))

# zap floats cache
_flt = float  # PYCHOK expected
del _Flts

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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

# % python -m pygeodesy.datum

# Ellipsoid.WGS84: name='WGS84', a=6378137, b=6356752.3142499998, f_=298.257223563, f=0.0033528107, e=0.0818191908, e2=0.00669438, e12=0.99330562, e22=0.0067394967, n=0.0016792204, R1=6371008.7714166669, R2=6371007.180920884, R3=6371000.7900107643, Rr=6367449.1458250275, Rs=6367435.6797186071,
#     A=6367449.145823414, e=8.1819190842621e-02, f_=1/298.2572235630, n=0.0016792203864(-3.791e-13),
#     AlphaKs=(0.00083773182062446983, 0.00000076085277735725, 0.00000000119764550324, 0.00000000000242917068, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#     BetaKs= (0.00083773216405794864, 0.0000000590587015222, 0.00000000016734826653, 0.00000000000021647981, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.00168275103155868244, 0.00504613293193328927, 0.00000529596776243447, 0.00000000690525779769)

# Ellipsoid.GRS80: name='GRS80', a=6378137, b=6356752.3141403468, f_=298.257222101, f=0.0033528107, e=0.081819191, e2=0.00669438, e12=0.99330562, e22=0.0067394968, n=0.0016792204, R1=6371008.7713801153, R2=6371007.1808835128, R3=6371000.7899741307, Rr=6367449.1457702452, Rs=6367435.6796636879,
#     A=6367449.145771047, e=8.1819191042816e-02, f_=1/298.2572221010, n=0.0016792203946(7.091e-16),
#     AlphaKs=(0.00083773182472855138, 0.00000076085278481497, 0.00000000119764552086, 0.00000000000242917073, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#     BetaKs= (0.00083773216816203529, 0.00000005905870210369, 0.000000000167348269, 0.00000000000021647982, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.00168275103983850371, 0.00504613295675166405, 0.0000052959678144849, 0.0000000069052578994)

# Ellipsoid.Sphere: name='Sphere', a=6371008.7714149999, b=6371008.7714149999, f_=0, f=0, e=0, e2=0, e12=1, e22=0, n=0, R1=6371008.7714149999, R2=6371008.7714149999, R3=6371008.7714149999, Rr=6371008.7714149999, Rs=6371008.7714149999,
#     A=6371008.771415, e=0.0e+00, f_=n/a, n=0.0(0.0e+00),
#     AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.0, 0.0, 0.0, 0.0)

# Ellipsoid.SpherePopular: name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, e=0, e2=0, e12=1, e22=0, n=0, R1=6378137, R2=6378137, R3=6378137, Rr=6378137, Rs=6378137,
#     A=6378137.0, e=0.0e+00, f_=n/a, n=0.0(0.0e+00),
#     AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.0, 0.0, 0.0, 0.0)
