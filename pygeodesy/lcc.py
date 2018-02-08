
# -*- coding: utf-8 -*-

u'''Lambert conformal conic projection for 1- or 2-Standard Parallels
class L{Conic}, L{Conics} registry and position class L{Lcc}.

See U{LCC<http://wikipedia.org/wiki/Lambert_conformal_conic_projection>},
U{Lambert Conformal Conic to Geographic Transformation Formulae
<http://www.linz.govt.nz/data/geodetic-system/coordinate-conversion/
projection-conversions/lambert-conformal-conic-geographic>},
U{Lambert Conformal Conic Projection
<http://mathworld.wolfram.com/LambertConformalConicProjection.html>}
and John P. Snyder U{'Map Projections - A Working Manual'
<http://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, 1987, pp 107-109.

@newfield example: Example, Examples
'''

from ellipsoidalBase import LatLonEllipsoidalBase as _eLLb
from datum import _Based, Datums, _Enum
from fmath import EPS, fStr, hypot
from utils import PI_2, degrees90, degrees180, false2f, \
                  radians, tanPI_2_2

from math import atan, copysign, cos, log, sin, sqrt

# all public constants, classes and functions
__all__ = ('Conic', 'Conics', 'Lcc',
           'toLcc')  # functions
__version__ = '18.02.06'


Conics = _Enum('Conics')  #: Registered conics (L{_Enum}).


class Conic(_Based):
    '''Lambert conformal conic projection (1- or 2-SP).
    '''
    _auth  = ''  #: (INTERNAL) authorization (string).
    _datum = None  #: (INTERNAL) Datum (L{Datum}).
    _name  = 'Conic'  #: (INTERNAL) Conic (L{Conic}).

    _e  = 0  #: (INTERNAL) Ellipsoid excentricity (float).
    _E0 = 0  #: (INTERNAL) False easting (float).
    _k0 = 1  #: (INTERNAL) Scale factor (float).
    _N0 = 0  #: (INTERNAL) false northing (float).
    _SP = 0  #: (INTERNAL) 1- or 2-SP (int)

    _lat0 = 0  #: (INTERNAL) Origin lat (radians).
    _lon0 = 0  #: (INTERNAL) Origin lon (radians).
    _par1 = 0  #: (INTERNAL) 1st std parallel (radians).
    _par2 = 0  #: (INTERNAL) 2nd std parallel (radians).
    _opt3 = 0  #: (INTERNAL) Optional, longitude (radians).

    _aF = 0  #: (INTERNAL) Precomputed F.
    _n  = 0  #: (INTERNAL) Precomputed n.
    _n_ = 0  #: (INTERNAL) Precomputed 1 / n.
    _r0 = 0  #: (INTERNAL) Precomputed rho0.

    def __init__(self, latlon0, par1, par2=None, E0=0, N0=0,
                       k0=1, opt3=0, name='', auth=''):
        '''New Lambert conformal conic projection.

           @param latlon0: Origin including an ellipsoidal datum (I{LatLon}).
           @param par1: First standard parallel (degrees90).
           @keyword par2: Optional, second standard parallel (degrees90).
           @keyword E0: Optional, false easting in meter (scalar).
           @keyword N0: Optional, false northing in meter (scalar).
           @keyword k0: Optional scale factor (scalar).
           @keyword opt3: Optional meridian (degrees180).
           @keyword name: Optional name of the conic (string).
           @keyword auth: Optional authentication authority (string).

           @return: A Lambert projection (L{Conic}).

           @raise TypeError: If I{latlon0} is not ellipsoidal.

           @example:

           >>> from pygeodesy import Conic, Datums, ellipsoidalNvector
           >>> ll0 = ellipsoidalNvector.LatLon(23, -96, datum=Datums.NAD27)
           >>> Snyder = Conic(ll0, 33, 45, E0=0, N0=0, name='Snyder')
        '''
        if latlon0 is not None:
            if not isinstance(latlon0, _eLLb):
                raise TypeError('%s not %s: %r' % ('latlon0', 'ellipsoidal', latlon0))

            self._lat0, self._lon0 = latlon0.to2ab()
            self._par1 = radians(par1)
            if par2 is None:
                self._par2 = self._par1
            else:
                self._par2 = radians(par2)
            self._opt3 = radians(opt3)

            if k0 != 1:
                self._k0 = float(k0)
            if N0:
                self._N0 = float(N0)
            if E0:
                self._E0 = float(E0)

            self.toDatum(latlon0.datum)._dup2(self)
            self._register(Conics, name)
        elif name:
            self._name = name
        if auth:
            self._auth = auth

    @property
    def auth(self):
        '''Get the authentication authority (string).
        '''
        return self._auth

    def copy(self, name=''):
        '''Copy this conic.

           @return: Conic, unregistered (L{Conic}).
        '''
        c = Conic(None, 0, name=name or self._name)
        self._dup2(c)
        return c

    @property
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property
    def E0(self):
        '''Get the false easting (meter).
        '''
        return self._E0

    @property
    def k0(self):
        '''Get scale factor (scalar).
        '''
        return self._k0

    @property
    def lat0(self):
        '''Get the origin latitude (degrees90).
        '''
        return degrees90(self._lat0)

    @property
    def lon0(self):
        '''Get the central meridian (degrees180).
        '''
        return degrees180(self._lon0)

    @property
    def N0(self):
        '''Get the false northing (meter).
        '''
        return self._N0

    @property
    def name(self):
        '''Get the conic name (string).
        '''
        return self._name

    @name.setter  # PYCHOK property setter
    def name(self, name):
        '''Set the conic name.

           @param name: New name (string).

           @raise ValueError: Invalid name.
        '''
        if name:
            self._name = str(name)
        else:
            raise ValueError('%s invalid: %r' % ('name', name))

    @property
    def name2(self):
        '''Get the conic and datum names as "conic.datum" (string).
        '''
        return self._name + '.' + self._datum.name

    @property
    def par1(self):
        '''Get the 1st standard parallel (degrees90).
        '''
        return degrees90(self._par1)

    @property
    def par2(self):
        '''Get the 2nd standard parallel (degrees90).
        '''
        return degrees90(self._par2)

    @property
    def opt3(self):
        '''Get the optional meridian (degrees180).
        '''
        return degrees180(self._opt3)

    @property
    def SP(self):
        '''Get the number of standard parallels (int).
        '''
        return self._SP

    def toDatum(self, datum):
        '''Convert this conic to the given datum.

           @param datum: Ellipsoidal datum to use (L{Datum}).

           @return: Converted conic, unregistered (L{Conic}).

           @raise TypeError: If I{datum} is not ellipsoidal.
        '''
        E = datum.ellipsoid
        if not E.isEllipsoidal:
            raise TypeError('%s not %s: %r' % ('datum', 'ellipsoidal', datum))

        c = self
        if c._e != E.e or c._datum != datum:

            c = self.copy()
            c._datum = datum

            c._e = E.e

            if abs(c._par1 - c._par2) < EPS:
                m1 = c._mdef(c._lat0)
                t1 = c._tdef(c._lat0)
                t0 = t1
                k  = 1
                n  = sin(c._lat0)
                sp = 1
            else:
                m1 = c._mdef(c._par1)
                m2 = c._mdef(c._par2)
                t1 = c._tdef(c._par1)
                t2 = c._tdef(c._par2)
                t0 = c._tdef(c._lat0)
                k  = c._k0
                n  = (log(m1) - log(m2)) \
                   / (log(t1) - log(t2))
                sp = 2

            F = m1 / (n * pow(t1, n))

            c._aF = k * E.a * F
            c._n  = n
            c._n_ = 1 / n
            c._r0 = c._rdef(t0)
            c._SP = sp

        return c

    convertDatum = toDatum  # alternate name

    def toStr(self, prec=8):  # PYCHOK expected
        '''Return this conic as a string.

           @keyword prec: Optional number of decimals, unstripped (int).

           @return: Conic attributes (string).
        '''
        if self._SP == 1:
            return self._fStr(prec, 'lat0', 'lon0', 'par1',
                                    'E0', 'N0', 'k0', 'SP',
                                     datum='(%s)' % (self.datum),)
        else:
            return self._fStr(prec, 'lat0', 'lon0', 'par1', 'par2',
                                    'E0', 'N0', 'k0', 'SP',
                                     datum='(%s)' % (self.datum),)

    def _dup2(self, c):
        '''(INTERNAL) Copy this conic to c.
           @param c: Duplicate (L{Conic}).
        '''
        c._auth  = self._auth
        c._datum = self._datum

        c._e  = self._e
        c._E0 = self._E0
        c._k0 = self._k0
        c._N0 = self._N0
        c._SP = self._SP

        c._lat0 = self._lat0
        c._lon0 = self._lon0
        c._par1 = self._par1
        c._par2 = self._par2
        c._opt3 = self._opt3

        c._aF = self._aF
        c._n  = self._n
        c._n_ = self._n_
        c._r0 = self._r0

    def _mdef(self, lat):
        '''(INTERNAL) Compute m(lat).
        '''
        s = self._e * sin(lat)
        return cos(lat) / sqrt(1 - s**2)

    def _pdef(self, lat):
        '''(INTERNAL) Compute p(lat).
        '''
        s = self._e * sin(lat)
        return pow((1 - s) / (1 + s), self._e / 2)

    def _rdef(self, t):
        '''(INTERNAL) Compute r(t).
        '''
        return self._aF * pow(t, self._n)

    def _tdef(self, lat):
        '''(INTERNAL) Compute t(lat).
        '''
        return max(0, tanPI_2_2(-lat) / self._pdef(lat))

    def _xdef(self, t_x):
        '''(INTERNAL) Compute x(t_x).
        '''
        return PI_2 - 2 * atan(t_x)  # XXX + self._lat0


Conics._assert(  # <http://spatialreference.org/ref/sr-org/...>
#   AsLb   = Conic(_eLLb(-14.2666667, 170, datum=Datums.NAD27), 0, 0, E0=500000, N0=0, name='AsLb', auth='EPSG:2155'),  # American Samoa ... SP=1 !
    Be08Lb = Conic(_eLLb(50.7978150, 4.359215833, datum=Datums.GRS80), 49.833333, 51.166667, E0=649328.0, N0=665262.0, name='Be08Lb', auth='EPSG:9802'),  # Belgium
    Be72Lb = Conic(_eLLb(90, 4.3674867, datum=Datums.NAD83), 49.8333339, 51.1666672, E0=150000.013, N0=5400088.438, name='Be72Lb', auth='EPSG:31370'),  # Belgium
    Fr93Lb = Conic(_eLLb(46.5, 3, datum=Datums.WGS84), 49, 44, E0=700000, N0=6600000, name='Fr93Lb', auth='EPSG:2154'),  # RFG93, France
    MaNLb  = Conic(_eLLb(33.3, -5.4, datum=Datums.NTF), 31.73, 34.87, E0=500000, N0=300000, name='MaNLb'),  # Marocco
    MxLb   = Conic(_eLLb(12, -102, datum=Datums.WGS84), 17.5, 29.5, E0=2500000, N0=0, name='MxLb', auth='EPSG:2155'),  # Mexico
    PyT_Lb = Conic(_eLLb(46.8, 2.33722917, datum=Datums.NTF), 45.89893890000052, 47.69601440000037, E0=600000, N0=200000, name='PyT_Lb', auth='Test'),  # France?
    USA_Lb = Conic(_eLLb(23, -96, datum=Datums.WGS84), 33, 45, E0=0, N0=0, name='USA_Lb'),  # Conterminous, contiguous USA?
    WRF_Lb = Conic(_eLLb(40, -97, datum=Datums.WGS84), 33, 45, E0=0, N0=0, name='WRF_Lb', auth='EPSG:4326')  # World
)


class Lcc(_Based):
    '''Lambert conformal conic East-/Northing location.
    '''
    _easting  = 0  #: (INTERNAL) Easting (float).
    _height   = 0  #: (INTERNAL) Height (meter).
    _northing = 0  #: (INTERNAL) Northing (float).
    _conic = None  #: (INTERNAL) Lamber projection (L{Conic}).

    def __init__(self, e, n, h=0, conic=Conics.WRF_Lb):
        '''New L{Lcc} position.

           @param e: Easting in meter (scalar).
           @param n: Northing in meter (scalar).
           @keyword h: Optional height in meter (scalar).
           @keyword conic: Optional, the conic projection (L{Conic}).

           @return: The Lambert location (L{Lcc}).

           @raise TypeError: If I{conic} is not L{Conic}.

           @raise ValueError: If I{e} or I{n} is invalid or negative.

           @example:

           >>> lb = Lcc(448251, 5411932.0001)
        '''
        if not isinstance(conic, Conic):
            raise TypeError('%s not Conic: %r' % ('conic', conic))
        self._conic = conic
        self._easting  = false2f(e, 'easting',  false=conic.E0 > 0)
        self._northing = false2f(n, 'northing', false=conic.N0 > 0)
        if h:
            self._height = float(h)

    @property
    def conic(self):
        '''Get the conic projection (L{Conic}).
        '''
        return self._conic

    @property
    def easting(self):
        '''Get the easting (meter).
        '''
        return self._easting

    @property
    def height(self):
        '''Get the height (meter).
        '''
        return self._height

    @property
    def northing(self):
        '''Get the northing (meter).
        '''
        return self._northing

    def to3lld(self, datum=None):
        '''Convert this L{Lcc} to a geodetic lat- and longitude.

           @keyword datum: Optional datum to use, otherwise use this
                           I{Lcc}'s conic.datum (I{Datum}).

           @return: 3-Tuple (lat, lon, datum) in (degrees90,
                    degrees180, datum).

           @raise TypeError: If I{datum} is not ellipsoidal.
        '''
        c = self.conic
        if datum:
            c = c.toDatum(datum)

        e =         self.easting  - c._E0
        n = c._r0 - self.northing + c._N0

        r_ = copysign(hypot(e, n), c._n)
        t_ = pow(r_ / c._aF, c._n_)

        x = c._xdef(t_)  # XXX c._lon0
        while True:
            p, x = x, c._xdef(t_ * c._pdef(x))
            if abs(x - p) < 1e-9:  # XXX EPS too small?
                break
        y = (atan(e / n) + c._opt3) * c._n_ + c._lon0

        return degrees90(x), degrees180(y), c.datum

    def toLatLon(self, LatLon, datum=None, height=None):
        '''Convert this L{Lcc} to an (ellipsoidal) geodetic point.

           @param LatLon: Ellipsoidal LatLon class to use for the
                          geodetic point (I{LatLon}).
           @keyword datum: Optional datum to use, otherwise use this
                           I{Lcc}'s conic.datum (I{Datum}).
           @keyword height: Optional height for the point, overriding
                            the default height (meter).

           @return: The point (I{LatLon}).

           @raise TypeError: If I{LatLon} or I{datum} is not ellipsoidal.
        '''
        if not issubclass(LatLon, _eLLb):
            raise TypeError('%s not %s: %r' % ('LatLon', 'ellipsoidal', LatLon))

        a, b, d = self.to3lld(datum=datum)

        h = self.height if height is None else height
        return LatLon(a, b, height=h, datum=d)

    def toStr(self, prec=0, sep=' ', m='m'):  # PYCHOK expected
        '''Return a string representation of this L{Lcc} position.

           @keyword prec: Optional number of decimal, unstripped (int).
           @keyword sep: Optional separator to join (string).
           @keyword m: Optional height units, default meter (string).

           @return: This Lcc as string "easting nothing" in meter plus
                    " height" and 'm' if heigth is non-zero (string).

           @example:

           >>> lb = Lcc(448251, 5411932.0001)
           >>> lb.toStr(4)  # 448251.0 5411932.0001
           >>> lb.toStr(sep=', ')  # 448251, 5411932
        '''
        t = [fStr(self._easting, prec=prec),
             fStr(self._northing, prec=prec)]
        if self._height:
            t += ['%+.2f%s' % (self._height, m)]
        return sep.join(t)

    def toStr2(self, prec=0, fmt='[%s]', sep=', ', m='m', C=False):  # PYCHOK expected
        '''Return a string representation of this L{Lcc} position.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword fmt: Optional, enclosing backets format (string).
           @keyword sep: Optional separator between name:values (string).
           @keyword m: Optional unit of the height, default meter (string).
           @keyword C: Optionally, include name of conic and datum (bool).

           @return: This Lcc as "[E:meter, N:meter, H:m, C:Conic.Datum]"
                   (string).
        '''
        t = self.toStr(prec=prec, sep=' ', m=m).split()
        k = 'ENH'[:len(t)]
        if C:
            k += 'C'
            t += [self.conic.name2]
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)


def toLcc(latlon, conic=Conics.WRF_Lb, height=None, Lcc=Lcc):
    '''Convert an (ellipsoidal) geodetic point to a Lambert location.

       @param latlon: Ellipsoidal point (I{LatLon}).
       @keyword conic: Optional Lambert projection to use (L{Conic}).
       @keyword height: Optional height for the point, overriding
                        the default height (meter).
       @keyword Lcc: Optional Lcc class to use for the Lambert
                     location (L{Lcc}).

       @return: The Lambert location (L{Lcc}).

       @raise TypeError: If I{latlon} is not ellipsoidal.
    '''
    if not isinstance(latlon, _eLLb):
        raise TypeError('%s not %s: %r' % ('latlon', 'ellipsoidal', latlon))

    c = conic.toDatum(latlon.datum)

    lat, lon = latlon.to2ab()
    r = c._rdef(c._tdef(lat))
    t = c._n * (lon - c._lon0) - c._opt3

    e = c._E0         + r * sin(t)
    n = c._N0 + c._r0 - r * cos(t)

    h = latlon.height if height is None else height
    return Lcc(e, n, h=h, conic=c)


if __name__ == '__main__':

    # print all
    for c in (Conics,):
        c = '\n' + repr(c)
        print('\n# '.join(c.split('\n')))

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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

# Conics.Be08Lb: Conic(name='Be08Lb', lat0=50.797815, lon0=4.35921583, par1=49.833333, par2=51.166667, E0=649328, N0=665262, k0=1, SP=2, datum=(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
# Conics.Be72Lb: Conic(name='Be72Lb', lat0=90, lon0=4.3674867, par1=49.8333339, par2=51.1666672, E0=150000.013, N0=5400088.438, k0=1, SP=2, datum=(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
# Conics.Fr93Lb: Conic(name='Fr93Lb', lat0=46.5, lon0=3, par1=49, par2=44, E0=700000, N0=6600000, k0=1, SP=2, datum=(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
# Conics.MaNLb: Conic(name='MaNLb', lat0=33.3, lon0=-5.4, par1=31.73, par2=34.87, E0=500000, N0=300000, k0=1, SP=2, datum=(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
# Conics.MxLb: Conic(name='MxLb', lat0=12, lon0=-102, par1=17.5, par2=29.5, E0=2500000, N0=0, k0=1, SP=2, datum=(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
# Conics.PyT_Lb: Conic(name='PyT_Lb', lat0=46.8, lon0=2.33722917, par1=45.8989389, par2=47.6960144, E0=600000, N0=200000, k0=1, SP=2, datum=(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
# Conics.USA_Lb: Conic(name='USA_Lb', lat0=23, lon0=-96, par1=33, par2=45, E0=0, N0=0, k0=1, SP=2, datum=(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
# Conics.WRF_Lb: Conic(name='WRF_Lb', lat0=40, lon0=-97, par1=33, par2=45, E0=0, N0=0, k0=1, SP=2, datum=(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)
