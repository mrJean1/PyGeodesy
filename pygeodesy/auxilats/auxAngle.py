
# -*- coding: utf-8 -*-

u'''(INTERNAL) I{Auxiliary} latitudes' base classes, constants and functions.

Class L{AuxAngle} transcoded to Python from I{Karney}'s C++ class U{AuxAngle
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxAngle.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2024) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.auxilats.auxily import Aux, _Aux2Greek,  AuxError
from pygeodesy.basics import map1, map2, _xinstanceof
from pygeodesy.constants import EPS, _INF_NAN_NINF, MAX, NAN, _0_0, _0_5, _1_0, \
                               _copysign_1_0, _over, _pos_self, isfinite, isnan
# from pygeodesy.errors import AuxError  # from .auxilats.auxily
from pygeodesy.fmath import hypot,  unstr
from pygeodesy.fsums import _add_op_, _iadd_op_, _isub_op_, _sub_op_
from pygeodesy.named import _Named,  _ALL_DOCS, _MODS
# from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS  # from .named
from pygeodesy.props import Property, Property_RO, property_RO, property_ROver, \
                           _update_all
# from pygeodesy.streprs import unstr  # from .fmath
from pygeodesy.units import Degrees, Radians
from pygeodesy.utily import atan2, atan2d, sincos2, sincos2d

from math import asinh, copysign, degrees, fabs, radians, sinh

__all__ = ()
__version__ = '25.05.12'

_0_INF_NAN_NINF = (0, _0_0) + _INF_NAN_NINF
_MAX_2          =  MAX * _0_5  # PYCHOK used!
# del _INF_NAN_NINF, MAX


class AuxAngle(_Named):
    '''U{An accurate representation of angles
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxAngle.html>}
    '''
    _AUX  =  None  # overloaded/-ridden
    _diff =  NAN   # default
    _iter =  None  # like ._NamedBase
    _y    = _0_0
    _x    = _1_0

    def __init__(self, y_angle=_0_0, x=_1_0, aux=None, **name):
        '''New L{AuxAngle}.

           @kwarg y_angle: The Y component (C{scalar}, including C{INF}, C{NAN}
                           and C{NINF}) or a previous L{AuxAngle} instance.
           @kwarg x: The X component, required if C{B{y_angle}} is C{scalar},
                     ignored otherwise.
           @kwarg aux: I{Auxiliary} kind (C{Aux.KIND}), like B{C{x}}.
           @kwarg name: Optional C{B{name}=NN} see (C{str}).

           @raise AuxError: Invalid B{C{y_angle}}, B{C{x}} or B{C{aux}}.
        '''
        try:
            try:
                yx  = y_angle._yx
                aux = y_angle._AUX
                if self._diff != y_angle._diff:
                    self._diff = y_angle._diff
            except AttributeError:
                yx  = y_angle, x
            if aux in _AUXClass:
                if self._AUX != aux:
                    self._AUX = aux
            elif aux is not None:
                raise ValueError()  # _invalid_
        except Exception as X:
            raise AuxError(y=y_angle, x=x, aux=aux, cause=X)
        self._y, self._x = _yx2(yx)
        if name:
            self.name = name

    def __abs__(self):
        '''Return this angle's absolute value (L{AuxAngle}).
        '''
        a     = self._copy_2(self.__abs__)
        a._yx = map2(fabs, self._yx)
        return a

    def __add__(self, other):
        '''Return C{B{self} + B{other}} as an L{AuxAngle}.

           @arg other: An L{AuxAngle}.

           @return: The sum (L{AuxAngle}).

           @raise TypeError: Invalid B{C{other}}.
        '''
        a = self._copy_2(self.__add__)
        return a._iadd(other, _add_op_)

    def __bool__(self):  # PYCHOK not special in Python 2-
        '''Return C{True} if this angle is non-zero.
        '''
        return bool(self.tan)

    def __eq__(self, other):
        '''Return C{B{self} == B{other}} as C{bool}.
        '''
        return not self.__ne__(other)

    def __float__(self):
        '''Return this angle's C{tan} (C{float}).
        '''
        return self.tan

    def __iadd__(self, other):
        '''Apply C{B{self} += B{other}} to this angle.

           @arg other: An L{AuxAngle}.

           @return: This angle, updated (L{AuxAngle}).

           @raise TypeError: Invalid B{C{other}}.
        '''
        return self._iadd(other, _iadd_op_)

    def __isub__(self, other):
        '''Apply C{B{self} -= B{other}} to this angle.

           @arg other: An L{AuxAngle}.

           @return: This instance, updated (L{AuxAngle}).

           @raise TypeError: Invalid B{C{other}} type.
        '''
        return self._iadd(-other, _isub_op_)

    def __ne__(self, other):
        '''Return C{B{self} != B{other}} as C{bool}.
        '''
        _xinstanceof(AuxAngle, other=other)
        y, x, r =  self._yxr_normalized()
        s, c, t = other._yxr_normalized()
        return fabs(y - s) > EPS or fabs(x - c) > EPS \
                                 or fabs(r - t) > EPS

    def __neg__(self):
        '''Return I{a copy of} this angle, negated.
        '''
        a = self._copy_2(self.__neg__)
        if a.y or not a.x:
            a.y = -a.y
        else:
            a.x = -a.x
        return a

    def __pos__(self):
        '''Return this angle I{as-is}, like C{float.__pos__()}.
        '''
        return self if _pos_self else self._copy_2(self.__pos__)

    def __radd__(self, other):
        '''Return C{B{other} + B{self}} as an L{AuxAngle}.

           @see: Method L{AuxAngle.__add__}.
        '''
        a = self._copy_r2(other, self.__radd__)
        return a._iadd(self, _add_op_)

    def __rsub__(self, other):
        '''Return C{B{other} - B{self}} as an L{AuxAngle}.

           @see: Method L{AuxAngle.__sub__}.
        '''
        a = self._copy_r2(other, self.__rsub__)
        return a._iadd(-self, _sub_op_)

    def __str__(self):
        n = _Aux2Greek.get(self._AUX, self.classname)
        return unstr(n, y=self.y, x=self.x, tan=self.tan)

    def __sub__(self, other):
        '''Return C{B{self} - B{other}} as an L{AuxAngle}.

           @arg other: An L{AuxAngle}.

           @return: The difference (L{AuxAngle}).

           @raise TypeError: Invalid B{C{other}} type.
        '''
        a = self._copy_2(self.__sub__)
        return a._iadd(-other, _sub_op_)

    def _iadd(self, other, *unused):  # op
        '''(INTERNAL) Apply C{B{self} += B{other}}.
        '''
        _xinstanceof(AuxAngle, other=other)
        # ignore zero other to preserve signs of _y and _x
        if other.tan:
            s, c = other._yx
            y, x = self._yx
            self._yx = (y * c + x * s), \
                       (x * c - y * s)
        return self

    def _copy_2(self, which):
        '''(INTERNAL) Copy for I{dyadic} operators.
        '''
        return _Named.copy(self, deep=False, name__=which)

    def _copy_r2(self, other, which):
        '''(INTERNAL) Copy for I{reverse-dyadic} operators.
        '''
        _xinstanceof(AuxAngle, other=other)
        return other._copy_2(which)

    def copyquadrant(self, other):
        '''Copy an B{C{other}} angle's quadrant into this angle (L{auxAngle}).
        '''
        _xinstanceof(AuxAngle, other=other)
        self._yx = copysign(self.y, other.y), \
                   copysign(self.x, other.x)
        return self

    @Property_RO
    def diff(self):
        '''Get derivative C{dtan(Eta) / dtan(Phi)} (C{float} or C{NAN}).
        '''
        return self._diff

    @staticmethod
    def fromDegrees(deg, **aux_name):
        '''Get an L{AuxAngle} from degrees.
        '''
        return _AuxClass(**aux_name)(*sincos2d(deg), **aux_name)

    @staticmethod
    def fromLambertianDegrees(psi, **aux_name):
        '''Get an L{AuxAngle} from I{Lambertian} degrees.
        '''
        return _AuxClass(**aux_name)(sinh(radians(psi)), **aux_name)

    @staticmethod
    def fromLambertianRadians(psi, **aux_name):
        '''Get an L{AuxAngle} from I{Lambertian} radians.
        '''
        return _AuxClass(**aux_name)(sinh(psi), **aux_name)

    @staticmethod
    def fromRadians(rad, **aux_name):
        '''Get an L{AuxAngle} from radians.
        '''
        return _AuxClass(**aux_name)(*sincos2(rad), **aux_name)

    @Property_RO
    def iteration(self):
        '''Get the iteration (C{int} or C{None}).
        '''
        return self._iter

    def normal(self):
        '''Normalize this angle I{in-place}.

           @return: This angle, normalized (L{AuxAngle}).
        '''
        self._yx = self._yx_normalized
        return self

    @Property_RO
    def normalized(self):
        '''Get a normalized copy of this angle (L{AuxAngle}).
        '''
        y, x = self._yx_normalized
        return self.classof(y, x, name=self.name, aux=self._AUX)

    @property_ROver
    def _RhumbAux(self):
        '''(INTERNAL) Import the L{RhumbAux} class, I{once}.
        '''
        return _MODS.rhumb.aux_.RhumbAux  # overwrite property_ROver

    @Property_RO
    def tan(self):
        '''Get this angle's C{tan} (C{float}).
        '''
        y, x = self._yx
        return _over(y, x) if isfinite(y) and y else y

    def toBeta(self, rhumb):
        '''Short for C{rhumb.auxDLat.convert(Aux.BETA, self, exact=rhumb.exact)}
        '''
        return self._toRhumbAux(rhumb, Aux.BETA)

    def toChi(self, rhumb):
        '''Short for C{rhumb.auxDLat.convert(Aux.CHI, self, exact=rhumb.exact)}
        '''
        return self._toRhumbAux(rhumb, Aux.CHI)

    @Property_RO
    def toDegrees(self):
        '''Get this angle as L{Degrees}.
        '''
        return Degrees(atan2d(*self._yx), name=self.name)

    @Property_RO
    def toLambertianDegrees(self):  # PYCHOK no cover
        '''Get this angle's I{Lambertian} in L{Degrees}.
        '''
        r = self.toLambertianRadians
        return Degrees(degrees(r), name=r.name)

    @Property_RO
    def toLambertianRadians(self):
        '''Get this angle's I{Lambertian} in L{Radians}.
        '''
        return Radians(asinh(self.tan), name=self.name)

    def toMu(self, rhumb):
        '''Short for C{rhumb.auxDLat.convert(Aux.MU, self, exact=rhumb.exact)}
        '''
        return self._toRhumbAux(rhumb, Aux.MU)

    def toPhi(self, rhumb):
        '''Short for C{rhumb.auxDLat.convert(Aux.PHI, self, exact=rhumb.exact)}
        '''
        return self._toRhumbAux(rhumb, Aux.PHI)

    @Property_RO
    def toRadians(self):
        '''Get this angle as L{Radians}.
        '''
        return Radians(atan2(*self._yx), name=self.name)

    def _toRhumbAux(self, rhumb, aux):
        '''(INTERNAL) Create an C{aux}-KIND angle from this angle.
        '''
        _xinstanceof(self._RhumbAux, rhumb=rhumb)
        return rhumb._auxD.convert(aux, self, exact=rhumb.exact)

    @Property
    def x(self):
        '''Get this angle's C{x} (C{float}).
        '''
        return self._x

    @x.setter  # PYCHOK setter!
    def x(self, x):  # PYCHOK no cover
        '''Set this angle's C{x} (C{float}).
        '''
        x = float(x)
        if self._x != x:
            _update_all(self)
            self._x = x

    @property_RO
    def _x_normalized(self):
        '''(INTERNAL) Get the I{normalized} C{x}.
        '''
        _, x = self._yx_normalized
        return x

    @Property
    def y(self):
        '''Get this angle's C{y} (C{float}).
        '''
        return self._y

    @y.setter  # PYCHOK setter!
    def y(self, y):  # PYCHOK no cover
        '''Set this angle's C{y} (C{float}).
        '''
        y = float(y)
        if self.y != y:
            _update_all(self)
            self._y = y

    @Property
    def _yx(self):
        '''(INTERNAL) Get this angle as 2-tuple C{(y, x)}.
        '''
        return self._y, self._x

    @_yx.setter  # PYCHOK setter!
    def _yx(self, yx):
        '''(INTERNAL) Set this angle's C{y} and C{x}.
        '''
        yx = _yx2(yx)
        if self._yx != yx:
            _update_all(self)
            self._y, self._x = yx

    @Property_RO
    def _yx_normalized(self):
        '''(INTERNAL) Get this angle as 2-tuple C{(y, x)}, I{normalized}.
        '''
        y, x = self._yx
        if isfinite(y) and fabs(y) < _MAX_2 \
                       and fabs(x) < _MAX_2 \
                       and isfinite(self.tan):
            h = hypot(y, x)
            if h > 0 and y:
                y = y / h  # /= chokes PyChecker
                x = x / h
                if isnan(y):  # PYCHOK no cover
                    y = _copysign_1_0(self.y)
                if isnan(x):  # PYCHOK no cover
                    x = _copysign_1_0(self.x)
            else:  # scalar 0
                y, x = _0_0, _copysign_1_0(y * x)
        else:  # scalar NAN
            y, x = NAN, _copysign_1_0(y * x)
        return y, x

    def _yxr_normalized(self, abs_y=False):
        '''(INTERNAL) Get 3-tuple C{(y, x, r)}, I{normalized}
           with C{y} or C{abs(y)} and C{r} as C{.toRadians}.
        '''
        y, x = self._yx_normalized
        if abs_y:
            y = fabs(y)  # only y, not x
        return y, x, atan2(y, x)  # .toRadians


class AuxBeta(AuxAngle):
    '''A I{Parametric, Auxiliary} latitude.
    '''
    _AUX = Aux.BETA

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxBeta} from degrees.
        '''
        return AuxBeta(*sincos2d(deg), **name)

    @staticmethod
    def fromRadians(rad, **name):
        '''Get an L{AuxBeta} from radians.
        '''
        return AuxBeta(*sincos2(rad), **name)


class AuxChi(AuxAngle):
    '''A I{Conformal, Auxiliary} latitude.
    '''
    _AUX = Aux.CHI

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxChi} from degrees.
        '''
        return AuxChi(*sincos2d(deg), **name)


class AuxMu(AuxAngle):
    '''A I{Rectifying, Auxiliary} latitude.
    '''
    _AUX = Aux.MU

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxMu} from degrees.
        '''
        return AuxMu(*sincos2d(deg), **name)


class AuxPhi(AuxAngle):
    '''A I{Geodetic or Geographic, Auxiliary} latitude.
    '''
    _AUX  =  Aux.PHI
    _diff = _1_0  # see .auxLat._Newton

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxPhi} from degrees.
        '''
        return AuxPhi(*sincos2d(deg), **name)


class AuxTheta(AuxAngle):
    '''A I{Geocentric, Auxiliary} latitude.
    '''
    _AUX = Aux.THETA

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxTheta} from degrees.
        '''
        return AuxTheta(*sincos2d(deg), **name)


class AuxXi(AuxAngle):
    '''An I{Authalic, Auxiliary} latitude.
    '''
    _AUX = Aux.XI

    @staticmethod
    def fromDegrees(deg, **name):
        '''Get an L{AuxXi} from degrees.
        '''
        return AuxXi(*sincos2d(deg), **name)


_AUXClass = {Aux.BETA:  AuxBeta,
             Aux.CHI:   AuxChi,
             Aux.MU:    AuxMu,
             Aux.PHI:   AuxPhi,
             Aux.THETA: AuxTheta,
             Aux.XI:    AuxXi}

def _AuxClass(aux=None, **unused):  # PYCHOK C{classof(aux)}
    return _AUXClass.get(aux, AuxAngle)


def _yx2(yx):
    '''(INTERNAL) Validate 2-tuple C{(y, x)}.
    '''
    try:
        y, x = yx
        y, x = map1(float, y, x)
        if y in _0_INF_NAN_NINF:
            x = _copysign_1_0(x)
    except (TypeError, ValueError) as X:
        raise AuxError(y=y, x=x, cause=X)
    return y, x


__all__ += _ALL_DOCS(AuxAngle, AuxBeta, AuxChi, AuxMu, AuxPhi, AuxTheta, AuxXi)

# **) MIT License
#
# Copyright (C) 2023-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
