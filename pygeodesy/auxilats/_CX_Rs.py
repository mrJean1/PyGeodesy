
# -*- coding: utf-8 -*-

u'''(INTERNAL) Classes C{_Rcoeffs}, C{_Rdict} and C{_Rtuple} to store the deferred
Python versions of coefficients from I{Karney}'s C++ class U{AuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}.

Copyright (C) Charles Karney (2022-2024) Karney@Alum.MIT.edu> and licensed under the
MIT/X11 License.  For more information, see <https:#GeographicLib.SourceForge.io>.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.constants import _floats
from pygeodesy.errors import _AssertionError,  _MODS
from pygeodesy.interns import NN, MISSING, _COMMA_, _duplicate_, _NL_, \
                             _QUOTE3_, _SLASH_,  _ELLIPSIS4_  # PYCHOK used!
# from pygeodesy.lazily import _ALL_MODS as _MODS  # from .errors
from pygeodesy.named import ADict,  Property_RO
# from pygeodesy.props import Property_RO  # from .named

__all__ = ()
__version__ = '24.09.04'


class _Rcoeffs(ADict):
    '''(INTERNAL) With string-ified C{keys}.
    '''
    def __init__(self, ALorder, coeffs):
        '''New C{_Rcoeffs} from a C{coeffs} dict.
        '''
        try:
            if not isinstance(coeffs, dict):
                raise _RdictError(coeffs=type(coeffs))
            n = 0
            for k, d in coeffs.items():
                if not isinstance(d, _Rdict):
                    raise _RdictError(k=k, d=type(d))
                n += d.n

            ADict.__init__(self, coeffs)
            self.set_(ALorder=ALorder, n=n)  # in .validate
        except Exception as x:
            raise _RdictError(ALorder=ALorder, cause=x)

    def bnuz4(self):  # in .auxilats.__main__  # PYCHOK no cover
        # get C{(strB, number, unique, floatB)} rationals
        b = n = u = z = 0
        _zB = _MODS.internals._sizeof
        for R in self._Rtuples():
            _, _, rs = R.k_n_rs
            b += _zB(rs)
            t  =  R._tuple
            z += _zB(t)  # Float
            # assert R.Rdict is None
            n += len(t)
            u += sum(1 for f in t if f in _floats)
        return b, n, (n - u), z

    def items(self):  # string-ify keys  # PYCHOK no cover
        for n, v in ADict.items(self):
            yield str(n), v

    def _Rtuples(self):  # PYCHOK no cover
        for d in self.values():
            if isinstance(d, _Rdict):
                # yield from d.values()
                for R in d.values():
                    yield R

    def _validate(self, aL, lenAux):
        # in .auxily.Aux._CXcoeffs(al, Aux.len(aL))
        a, n = self.ALorder, self.n  # PYCHOK Adict!
#       for R in self._Rtuples():
#           assert isinstance(R, _Rtuple)
        if aL != a or lenAux != n:
            raise _RdictError(aL=aL, ALorder=a, lenAux=lenAux, n=n)
        return self


class _Rdict(dict):  # in ._CX_#, .auxLat, .rhumb.aux_
    '''(INTERNAL) Dict of C{_Rtuple}s.
    '''
    n = 0  # sum(R.k_n_k[1] for r in Rtuples)

    def __init__(self, nt, *Rtuples):
        '''New C{_Rdict}.
        '''
        if not Rtuples:
            raise _RdictError(Rtuples=MISSING)

        for R in Rtuples:
            if not isinstance(R, _Rtuple):
                raise _RdictError(R, R=type(R))
            k, n, _ = R.k_n_rs
            if k in self:
                raise _RdictError(_duplicate_, k=k)
            R.Rdict = self
            self[k] = R  # overwritten in self._floatuple
            self.n += n
        if self.n != nt:
            raise _RdictError(n=n, nt=nt)

    def _floats(self, rs):
        # yield floats from a string of comma-separated rationals
        def _p_q(p=NN, q=1, *x):
            return (NN if x else p), q

        _get = _floats.get
        for r in NN(*rs.split()).split(_COMMA_):
            p, q = _p_q(*r.split(_SLASH_))
            if p:
                f = int(p) / int(q)  # fractions.Fraction?
                if not isinstance(f, float):
                    raise _RdictError(rs, f=f, p=p, q=q, r=r)
                yield _get(f, f)  # from .constants?
            else:
                raise _RdictError(rs, r=r)

    def _floatuple(self, Rtuple):
        # return a tuple of floats from an C{_Rtuple}
        k, n, rs = Rtuple.k_n_rs
        t = tuple(f for m in map(self._floats, rs)
                        for f in m)  # ... yield f
        # @see: <https://StackOverflow.com/questions/10632839/>
        #       and <https://realPython.com/python-flatten-list/>
        if len(t) != n:
            raise _RdictError(*rs, len=len(t), n=n)
        self[k] = t
        return t


class _RdictError(_AssertionError):
    '''(INTERNAL) For C{_Rdict} issues.
    '''
    def __init__(self, *rs, **kwds_cause):  # PYCHOK no cover
        if rs:
            if len(rs) > 1:
                t = _NL_(NN, *rs)
                t =  NN(_QUOTE3_, t, _QUOTE3_)
            else:  # single rs
                t =  repr(rs[0])
            kwds_cause.update(txt=t)
        _AssertionError.__init__(self, **kwds_cause)


class _Rtuple(list):  # MUST be list, NOT tuple!
    '''(INTERNAL) I{Pseudo-tuple} of float rationals used in C{_Rdict}s.
    '''
    Rdict  = None
    k_n_rs = None, 0, ()

    def __init__(self, k, n, *rs):
        '''New C{_Rtuple} with key C{k}, number of floats C{n} and with
           each C{rs} a C{str} of comma-separated rationals C{"p/q, ..."}
           where C{p} and C{q} are C{int} digits only.
        '''
        try:
            if not rs:
                raise _RdictError(rs=MISSING)
            for t in rs:
                if not isinstance(t, str):
                    raise _RdictError(rs=type(t))
            if not (isinstance(n, int) and n > 0):
                raise _RdictError(n=type(n))
            self.k_n_rs = k, n, rs
        except Exception as x:
            raise _RdictError(*rs, k=k, n=n, cause=x)

    def __getitem__(self, i):
        return self._tuple[i]

    def __iter__(self):
        return iter(self._tuple)

    def __len__(self):
        return len(self._tuple)

    @Property_RO
    def _tuple(self):
        # build the C{tuple} once, replace C{_Rdict}
        # item at C{key} with the C{tuple} and fill
        # this C{_Rlist} with the C{tuple} values
        # for the initial __getitem__ retrieval[s]
        try:
            k, n, rs = self.k_n_rs
            t = self.Rdict._floatuple(self)
            self[:] = t  # MUST copy into self!
        except Exception as x:
            if len(rs) > 1 and _QUOTE3_ in str(x):
                rs = rs[0], _ELLIPSIS4_
            raise _RdictError(k=k, n=n, rs=rs, cause=x)
        del self.Rdict, self.k_n_rs  # trash refs
        return t

    def append(self, arg):
        raise _RdictError(append=arg)

    def extend(self, arg):
        raise _RdictError(extend=arg)

# **) MIT License
#
# Copyright (C) 2024-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
