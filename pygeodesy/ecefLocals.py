
# -*- coding: utf-8 -*-

u'''(INTERNAL) ECEF to local coodinate conversions, separated from
module C{ecef} to defer importing the latter into C{CartesianBase}
and C{LatLonBase}.
'''

from pygeodesy.basics import _isin, _xsubclassof
# from pygeodesy.ecef import EcefKarney  # _MODS
from pygeodesy.errors import _xkwds_item2, _xkwds_pop2
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy import ltp as _ltp  # _MODS.into
# from pygeodesy import ltpTuples as _ltpTuples  # _MODS.into
from pygeodesy.named import _Named, notOverloaded
from pygeodesy.props import Property_RO, property_RO, property_ROver

__all__ = _ALL_LAZY.ecefLocals
__version__ = '25.04.28'

_ltp       = _MODS.into(ltp=__name__)
_ltpTuples = _MODS.into(ltpTuples=__name__)


class _EcefLocal(_Named):
    '''(INTERNAL) Base class for C{CartesianBase}, C{Ecef9Tuple} and
       C{LatLonBase}, providing ECEF to local coordinate conversions.
    '''

    @property_ROver
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney}), I{once}.
        '''
        return _MODS.ecef.EcefKarney

    @property_RO
    def _ecef9(self):
        '''I{Must be overloaded}.'''
        notOverloaded(self)

    @property_RO
    def _ecef9datum(self):
        try:
            d = self.datum  # PYCHOK C{CartesianBase}, ...
        except AttributeError:
            d = None
        return d or self._ecef9.datum

    @Property_RO
    def _ltp(self):
        '''(INTERNAL) Cache this instance' LTP (L{Ltp}).
        '''
        return _ltp.Ltp(self._ecef9, ecef=self.Ecef(self._ecef9datum), name=self.name)

    def _ltp_ecef2local(self, ltp, Xyz_kwds, _None=None, **Xyz):  # in C{Ecef9Tuple}
        '''(INTERNAL) Invoke C{_ltp._xLtp(ltp, self._ltp)._ecef2local}.
        '''
        C, _ = Xyz_ = _xkwds_pop2(Xyz_kwds, **Xyz)
        if C is not _None:  # validate C, see .toLocal
            n, X = _xkwds_item2(Xyz)
            if X is not C:
                _xsubclassof(X, **{n: C})
        ltp = _ltp._xLtp(ltp, self._ltp)
        return ltp._ecef2local(self._ecef9, *Xyz_)

    def toAer(self, ltp=None, **Aer_and_kwds):
        '''Convert this instance to I{local} I{Azimuth, Elevation, slant Range} (AER) components.

           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}), overriding this
                       instance' L{LTP<pygeodesy.ecefLocals._EcefLocal.toLtp>}.
           @kwarg Aer_and_kwds: Optional AER class C{B{Aer}=}L{Aer<pygeodesy.ltpTuples.Aer>}
                      to use and optionally, additional B{C{Aer}} keyword arguments.

           @return: An B{C{Aer}} instance.

           @raise TypeError: Invalid B{C{ltp}}.

           @see: Method L{toLocal<pygeodesy.ecefLocals._EcefLocal.toLocal>}.
        '''
        return self._ltp_ecef2local(ltp, Aer_and_kwds, Aer=_ltpTuples.Aer)

    def toEnu(self, ltp=None, **Enu_and_kwds):
        '''Convert this instance to I{local} I{East, North, Up} (ENU) components.

           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}), overriding this
                       instance' L{LTP<pygeodesy.ecefLocals._EcefLocal.toLtp>}.
           @kwarg Enu_and_kwds: Optional ENU class C{B{Enu}=}L{Enu<pygeodesy.ltpTuples.Enu>}
                      to use and optionally, additional B{C{Enu}} keyword arguments.

           @return: An B{C{Enu}} instance.

           @raise TypeError: Invalid B{C{ltp}}.

           @see: Method L{toLocal<pygeodesy.ecefLocals._EcefLocal.toLocal>}.
        '''
        return self._ltp_ecef2local(ltp, Enu_and_kwds, Enu=_ltpTuples.Enu)

    def toLocal(self, Xyz=None, ltp=None, **Xyz_kwds):
        '''Convert this instance to I{local} components in a I{local tangent plane} (LTP)

           @kwarg Xyz: Optional I{local} components class (L{XyzLocal}, L{Aer}, L{Enu},
                       L{Ned}) or C{None}.
           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}), overriding this
                       cartesian's L{LTP<pygeodesy.ecefLocals._EcefLocal.toLtp>}.
           @kwarg Xyz_kwds: Optionally, additional B{C{Xyz}} keyword arguments, ignored
                            if C{B{Xyz} is None}.

           @return: An B{C{Xyz}} instance or a L{Local9Tuple}C{(x, y, z, lat, lon,
                    height, ltp, ecef, M)} if C{B{Xyz} is None} (with C{M=None}).

           @raise TypeError: Invalid B{C{ltp}}.
        '''
        return self._ltp_ecef2local(ltp, Xyz_kwds, Xyz=Xyz, _None=Xyz)

    def toLtp(self, Ecef=None, **name):
        '''Return the I{local tangent plane} (LTP) for this instance.

           @kwarg Ecef: Optional ECEF I{class} (L{EcefKarney}, ... L{EcefYou}), overriding
                        this instance' L{Ecef<pygeodesy.ecefLocals._EcefLocal.Ecef>}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: An B{C{Ltp}} instance.
        '''
        if _isin(Ecef, None, self.Ecef) and not name:
            ltp = self._ltp
        else:  # like self._ltp
            ecef = (Ecef if Ecef else self.Ecef)(self._ecef9datum)
            ltp  = _ltp.Ltp(self._ecef9, ecef=ecef, name=self._name__(name))
        return ltp

    def toNed(self, ltp=None, **Ned_and_kwds):
        '''Convert this instance to I{local} I{North, East, Down} (NED) components.

           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}), overriding this
                       instance' L{LTP<pygeodesy.ecefLocals._EcefLocal.toLtp>}.
           @kwarg Ned_and_kwds: Optional NED class C{B{Ned}=}L{Ned<pygeodesy.ltpTuples.Ned>}
                      to use and optionally, additional B{C{Ned}} keyword arguments.

           @return: An B{C{Ned}} instance.

           @raise TypeError: Invalid B{C{ltp}}.

           @see: Method L{toLocal<pygeodesy.ecefLocals._EcefLocal.toLocal>}.
        '''
        return self._ltp_ecef2local(ltp, Ned_and_kwds, Ned=_ltpTuples.Ned)

    def toXyz(self, ltp=None, **Xyz_and_kwds):
        '''Convert this instance to I{local} I{X, Y, Z} (XYZ) components.

           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}), overriding this
                       instance' L{LTP<pygeodesy.ecefLocals._EcefLocal.toLtp>}.
           @kwarg Xyz_and_kwds: Optional XYZ class C{B{Xyz}=}L{Xyz<pygeodesy.ltpTuples.XyzLocal>}
                      to use and optionally, additional B{C{Xyz}} keyword arguments.

           @return: An B{C{Xyz}} instance.

           @raise TypeError: Invalid B{C{ltp}}.

           @see: Method L{toLocal<pygeodesy.ecefLocals._EcefLocal.toLocal>}.
        '''
        return self._ltp_ecef2local(ltp, Xyz_and_kwds, Xyz=_ltpTuples.XyzLocal)


__all__ += _ALL_DOCS(_EcefLocal)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
