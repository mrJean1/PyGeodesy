
# -*- coding: utf-8 -*-

# Lambert Conformal Conic projection for 1 or 2 Standard Parallels.

# <http://wikipedia.org/wiki/Lambert_conformal_conic_projection>
# <http://www.linz.govt.nz/data/geodetic-system/coordinate-conversion/
#         projection-conversions/lambert-conformal-conic-geographic>
# Snyder <https://pubs.er.usgs.gov/djvu/PP/PP_1395.pdf> pp 107-109
# <http://mathworld.wolfram.com/LambertConformalConicProjection.html>

from math import atan, copysign, cos, hypot, log, sin, sqrt, tan
from ellipsoidalBase import _LatLonHeightDatumBase as _LL
from datum import _Based, Datums, _Enum
from utils import EPS, PI_2, degrees90, degrees180, fStr, radians

# all public constants, classes and functions
__all__ = ('Conic', 'Conics', 'Lcc',
           'toLcc')  # functions
__version__ = '16.12.19'


Conics = _Enum('Conics')


class Conic(_Based):
    '''Lambert conformal conic projection (1- or 2-SP).
    '''
    _auth  = ''
    _datum = None
    _name  = 'Conic'

    _e  = 0  # ellipsoid excentricity
    _E0 = 0  # false easting
    _k0 = 1  # scale factor
    _N0 = 0  # false northing
    _SP = 0  # 1- or 2-SP

    _lat0 = 0  # origin lat in radians
    _lon0 = 0  # origin lon in radians
    _par1 = 0  # 1st std parallel in radians
    _par2 = 0  # 2nd std parallel in radians
    _opt3 = 0  # XXX optional, longitude in radians

    _aF = 0  # precomputed values
    _n  = 0
    _n_ = 0
    _r0 = 0

    def __init__(self, latlon0, par1, par2, E0=0, N0=0,
                       k0=1, opt3=0, name='', auth=''):
        '''Create a Lambert conformal conic projection.

           @param {LatLon} latlon - Origin in degrees including
                                    an ellipsoidal datum.
           @param {degrees90} par1 - First standard parallel.
           @param {degrees90} par2 - Second standard parallel.
           @param {meter} [E0=0] - False easting in meter.
           @param {meter} [N0=0] - False northing in meter.
           @param {number} [k0=1] - Scale factor.
           @param {degrees180} [opt3=0] - Optional, longitude.
           @param {string} [name=''] - Name of the conic.
           @param {string} [auth=''] - Authentication authority.

           @returns {Conic} The Lambert projection.

           @example
           from geodesy import Conic, Datums, ellipsoidalNvector
           LL = ellipsoidalNvector.LatLon
           Snyder = Conic(LL(23, -96, datum=Datums.NAD27),
                          33, 45, E0=0, N0=0, name='Snyder')
        '''
        if latlon0 is not None:
            if not isinstance(latlon0, _LL):
                raise TypeError('%s not ellipsoidal: %r' % ('latlon0', latlon0))

            self._lat0, self._lon0 = latlon0.toradians()
            self._par1 = radians(par1)
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
        '''Return the authentication.
        '''
        return self._auth

    def copy(self, name=''):
        '''Return a copy of this Conic, unregistered.
        '''
        c = Conic(None, 0, 0, name=name or self._name)
        self._dup2(c)
        return c

    @property
    def datum(self):
        '''Return the datum.
        '''
        return self._datum

    @property
    def E0(self):
        '''Return false easting in meter.
        '''
        return self._E0

    @property
    def k0(self):
        '''Return scale factor.
        '''
        return self._k0

    @property
    def lat0(self):
        '''Return latitude of origin in degrees90.
        '''
        return degrees90(self._lat0)

    @property
    def lon0(self):
        '''Return central meridian in degrees180.
        '''
        return degrees180(self._lon0)

    @property
    def N0(self):
        '''Return false northing in meter.
        '''
        return self._N0

    @property
    def name(self):
        '''Return the name.
        '''
        return self._name

    @name.setter  # PYCHOK property setter
    def name(self, name):
        '''Set the conic name.
        '''
        self._name = name

    @property
    def name2(self):
        '''Return the conic and datum name.
        '''
        return self._name + '.' + self._datum.name

    @property
    def par1(self):
        '''Return 1st standard parallel in degrees90.
        '''
        return degrees90(self._par1)

    @property
    def par2(self):
        '''Return 2nd standard parallel in degrees90.
        '''
        return degrees90(self._par2)

    @property
    def opt3(self):
        '''Return optional longitude in degrees180.
        '''
        return degrees180(self._opt3)

    @property
    def SP(self):
        '''Return 1- or 2-SP.
        '''
        return self._SP

    def toDatum(self, datum):
        '''Convert this conic to the given datum.
        '''
        E = datum.ellipsoid
        if not E.b < E.R < E.a:
            raise ValueError('%s not ellipsoidal: %r' % ('datum', datum))

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

           @param {int} [prec=8] - Number of decimals, unstripped.

           @returns {string} Conic string.
        '''
        return self._fStr(prec, 'lat0', 'lon0', 'par1', 'par2',
                                'E0', 'N0', 'k0', 'SP',
                                 datum='(%s)' % (self.datum),)

    def _dup2(self, c):
        # copy self to c
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

    def _mdef(self, lat):  # compute m(lat)
        s = self._e * sin(lat)
        return cos(lat) / sqrt(1 - s * s)

    def _pdef(self, lat):  # compute p(lat)
        s = self._e * sin(lat)
        return pow((1 - s) / (1 + s), self._e / 2)

    def _rdef(self, t):  # compute r(t)
        return self._aF * pow(t, self._n)

    def _tdef(self, lat):  # compute t(lat)
        return max(0, tan((PI_2 - lat) / 2) / self._pdef(lat))

    def _xdef(self, t_x):
        return PI_2 - 2 * atan(t_x)  # XXX + self._lat0


Conics._assert(  # <http://spatialreference.org/ref/sr-org/...>
#   AsLb   = Conic(_LL(-14.2666667, 170, datum=Datums.NAD27), 0, 0, E0=500000, N0=0, name='AsLb', auth='EPSG:2155'),  # American Samoa ... SP=1 !
    Be72Lb = Conic(_LL(90, 4.3674867, datum=Datums.NAD83), 49.8333339, 51.1666672, E0=150000.013, N0=5400088.438, name='Be72Lb', auth='EPSG:31370'),  # Belgium
    Fr93Lb = Conic(_LL(46.5, 3, datum=Datums.WGS84), 49, 44, E0=700000, N0=6600000, name='Fr93Lb', auth='EPSG:2154'),  # RFG93, France
    MaNLb  = Conic(_LL(33.3, -5.4, datum=Datums.NTF), 31.73, 34.87, E0=500000, N0=300000, name='MaNLb'),  # Marocco
    MxLb   = Conic(_LL(12, -102, datum=Datums.WGS84), 17.5, 29.5, E0=2500000, N0=0, name='MxLb', auth='EPSG:2155'),  # Mexico
    PyT_Lb = Conic(_LL(46.8, 2.33722917, datum=Datums.NTF), 45.89893890000052, 47.69601440000037, E0=600000, N0=200000, name='PyT_Lb', auth='Test'),  # France?
    USA_Lb = Conic(_LL(23, -96, datum=Datums.WGS84), 33, 45, E0=0, N0=0, name='USA_Lb'),  # Conterminous, contiguous USA?
    WRF_Lb = Conic(_LL(40, -97, datum=Datums.WGS84), 33, 45, E0=0, N0=0, name='WRF_Lb', auth='EPSG:4326')  # World
)


class Lcc(_Based):
    '''Lambert conformal conic East-/Northing location.
    '''

    def __init__(self, e, n, h=0, conic=Conics.WRF_Lb):
        '''Create an Lcc position.

           @param {meter} e - Easting.
           @param {meter} n - Northing.
           @param {meter} [h=0] - Height
           @param {Conic} [conic=Conics.WRF_Lb] - The conic projection.

           @returns {Lcc} The Lcc location.

           @example
           lb = Lcc(448251, 5411932.0001)
        '''
        self._easting = e
        self._northing = n
        self._height = h
        self._conic = conic

    @property
    def conic(self):
        '''Return the conic projection.
        '''
        return self._conic

    @property
    def easting(self):
        '''Return the easting in meter.
        '''
        return self._easting

    @property
    def height(self):
        '''Return the height in meter.
        '''
        return self._height

    @property
    def northing(self):
        '''Return the northing in meter.
        '''
        return self._northing

    def toLatLon(self, LatLon=_LL, datum=None):
        '''Convert this Lcc to a lat-/longitude instance.

           @param {LatLon} LatLon - LatLon (sub-)class to instantiate.
           @param (Datum} [datum=None] - LatLon datum to use, otherwise
                                         use this Lcc's conic.datum.

           @returns {LatLon} The position as a LatLon instance.
        '''
        if not issubclass(LatLon, _LL):
            raise TypeError('%s not ellipsoidal: %r' % ('LatLon', LatLon))

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
            if abs(x - p) < 1e-9:  # EPS
                break
        y = (atan(e / n) + c._opt3) * c._n_ + c._lon0

        return LatLon(degrees90(x), degrees180(y), height=self.height, datum=c.datum)

    def toStr(self, prec=0, sep=' ', m='m'):  # PYCHOK expected
        '''Returns a string representation of this Lcc position.

           @param {number} [prec=0] - Number of decimal, unstripped.
           @param {string} [sep=' '] - Separator to join.
           @param {string} [m='m'] - Unit of the height, default meter.

           @returns {string} Lcc as string "easting nothing" in meter
                             plus " height" in m if non-zero.

           @example
           lb = Lcc(448251, 5411932.0001)
           lb.toStr(4)  # 448251.0 5411932.0001
           lb.toStr(sep=', ')  # 448251, 5411932
        '''
        t = [fStr(self._easting, prec=prec),
             fStr(self._northing, prec=prec)]
        if self._height:
            t += ['%+.2f%s' % (self._height, m)]
        return sep.join(t)

    def toStr2(self, prec=0, fmt='[%s]', sep=', ', m='m', C=False):  # PYCHOK expected
        '''Returns a string representation of this Lcc position.

           @param {number} [prec=0] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between name:values.
           @param {string} [m='m'] - Unit of the height, default meter.
           @param {bool} [C=False] - Include name of conic and datum.

           @returns {string} This Lcc as "[E:meter, N:meter, H:m, C:Conic.Datum]".
        '''
        t = self.toStr(prec=prec, sep=' ', m=m).split()
        k = 'ENH'[:len(t)]
        if C:
            k += 'C'
            t += [self.conic.name2]
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)


def toLcc(latlon, conic=Conics.WRF_Lb):
    '''Convert a lat-/longitude position to an Lcc instance.

       @param {LatLon} latlon - Lat-/longitude location.
       @parma {Conic} [conic=Conics.WRF_Lb] - Lambert conic to use.

       @returns {Lcc} The location as an Lcc instance.
    '''
    if not isinstance(latlon, _LL):
        raise TypeError('%s not ellipsoidal: %r' % ('latlon', latlon))

    c = conic.toDatum(latlon.datum)

    lat, lon = latlon.toradians()
    r = c._rdef(c._tdef(lat))
    t = c._n * (lon - c._lon0) - c._opt3

    return Lcc(c._E0         + r * sin(t),
               c._N0 + c._r0 - r * cos(t), h=latlon.height, conic=c)


if __name__ == '__main__':

    # print all
    for c in (Conics,):
        print('\n%r' % (c,))

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1@Gmail.com
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

# Typical result (on MacOS X)

# Conics.Be72Lb: Conic(lat0=90.0, lon0=4.3674867, par1=49.8333339, par2=51.1666672, E0=150000.013, N0=5400088.438, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83, name='NAD83'), name='Be72Lb')
# Conics.Fr93Lb: Conic(lat0=46.5, lon0=3.0, par1=49.0, par2=44.0, E0=700000.0, N0=6600000.0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84, name='WGS84'), name='Fr93Lb')
# Conics.MaNLb: Conic(lat0=33.3, lon0=-5.4, par1=31.73, par2=34.87, E0=500000.0, N0=300000.0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF, name='NTF'), name='MaNLb')
# Conics.MxLb: Conic(lat0=12.0, lon0=-102.0, par1=17.5, par2=29.5, E0=2500000.0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84, name='WGS84'), name='MxLb')
# Conics.PyT_Lb: Conic(lat0=46.8, lon0=2.33722917, par1=45.8989389, par2=47.6960144, E0=600000.0, N0=200000.0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF, name='NTF'), name='PyT_Lb')
# Conics.USA_Lb: Conic(lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84, name='WGS84'), name='USA_Lb')
# Conics.WRF_Lb: Conic(lat0=40.0, lon0=-97.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84, name='WGS84'), name='WRF_Lb')
