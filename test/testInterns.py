
# -*- coding: utf-8 -*-

# Test L{interns} module.

__all__ = ('Tests',)
__version__ = '23.08.22'

from bases import ismacOS, sys, TestsBase

from pygeodesy import clips, interns, machine, NN

from os import getcwd
# import sys  # from .bases

_0to9_      =  interns._0to9_
_AtoZnoIO_  =  interns._AtoZnoIO_
_cwd        =  getcwd()
_functions  = (interns._platform2,
               interns._version2)
_exceptions = (_0to9_, _AtoZnoIO_,
               interns._doesn_t_exist_,
               interns._exceed_PI_radians_,
               interns._iadd_op_,
               interns._n_a_,
               interns._NLATvar_,
               interns._NLHASH_,
               interns._not_finite_,
               interns._not_scalar_,
               interns._PyPy__,
               interns._semi_circular_,
               interns._utf_8_)


class Tests(TestsBase):

    def testInterns(self):

        _COLONSPACE_  = interns._COLONSPACE_
        _DUNDER_      = interns._DUNDER_
        _DOT_         = interns._DOT_
        _DUNDER_0to9_ = NN(_DUNDER_, _0to9_)
        _EQUALSPACED_ = interns._EQUALSPACED_
        _UNDER_       = interns._UNDER_

        for a in sorted(dir(interns), key=str.lower):
            if a.startswith(_UNDER_) and a[-1:] in _DUNDER_0to9_:
                i = getattr(interns, a, None)
                n = repr(i)
                k = a.startswith(_DUNDER_)
                if k:  # hide home dir
                    n = n.replace(_cwd, _DOT_)
                n = _EQUALSPACED_(a, clips(n))
                self.test(n, isinstance(i, (float, str)) or i in _functions, True, known=k)
                # check the naming conventions
                if isinstance(i, str) and not k:
                    a = a.strip(_UNDER_)
                    self.test(n, i.lower(), a.lower(), known=a.isupper() or i in _exceptions)
                elif isinstance(i, float):
                    a = a.strip(_UNDER_).replace(_UNDER_, _DOT_)
                    self.test(n, i, a, known=n.startswith(_UNDER_))

        self.test('.tillC', _AtoZnoIO_.tillC, 'ABC')
        self.test('.fromX', _AtoZnoIO_.fromX, 'XYZ')

        self.test('.fromH.tillJ', _AtoZnoIO_.fromH.tillJ, 'HJ')
        self.test('.fromN.tillP', _AtoZnoIO_.fromN.tillP, 'NP')

        m = machine()
        self.test('machine', m, m, nl=1)
        if ismacOS:  # coverage
            i = interns._sysctl_uint('sysctl.proc_translated')
            self.test('sysctl', i, i)  # known=i in (0, 1)
            u, t = interns._usage('interns.py').split(_COLONSPACE_)
            self.test(u, t, t)
            t = interns._version2(sys.version)
            self.test('version', t, tuple(sys.version_info[:2]))


if __name__ == '__main__':

    t = Tests(__file__, __version__, interns)
    t.testInterns()
    t.results()
    t.exit()
