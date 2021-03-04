
# -*- coding: utf-8 -*-

u'''Test L{interns} module.
'''

__all__ = ('Tests',)
__version__ = '21.01.24'

from base import TestsBase

from pygeodesy import clips, interns, EPS, EPS1, NN

from os import getcwd

_0to9_        =  interns._0to9_
_AtoZnoIO_    =  interns._AtoZnoIO_
_cwd          =  getcwd()
_DUNDER_      =  interns._DUNDER_
_DOT_         =  interns._DOT_
_EPS4         =  interns._EPS4
_EPS__2       =  interns._EPS__2
_EPS__4       =  interns._EPS__4
_EQUALSPACED_ =  interns._EQUALSPACED_
_UNDER_       =  interns._UNDER_
_exceptions   = (_0to9_, _AtoZnoIO_,
                 interns._doesn_t_exist_,
                 interns._exceed_PI_radians_,
                 interns._n_a_,
                 interns._near_concentric_,
                 interns._NL_hash_,
                 interns._NL_var_,
                 interns._OKd_,
                 interns._semi_circular_,
                 interns._utf_8_)
_DUNDER_0to9_ =  NN(_DUNDER_, _0to9_)


class Tests(TestsBase):

    def testInterns(self):

        for a in sorted(dir(interns), key=str.lower):
            if a.startswith(_UNDER_) and a[-1:] in _DUNDER_0to9_:
                i = getattr(interns, a, None)
                n = repr(i)
                k = a.startswith(_DUNDER_)
                if k:  # hide home dir
                    n = n.replace(_cwd, _DOT_)
                n = _EQUALSPACED_(a, clips(n))
                self.test(n, isinstance(i, (float, str)), True, known=k)
                # check the naming conventions
                if isinstance(i, str) and not k:
                    a = a.strip(_UNDER_)
                    self.test(n, i.lower(), a.lower(), known=a.isupper() or i in _exceptions)
                elif isinstance(i, float) and i not in (_EPS4, _EPS__2, _EPS__4):
                    a = a.strip(_UNDER_).replace(_UNDER_, _DOT_)
                    self.test(n, i, a)

        self.test('EPS',    EPS > 0, True)
        self.test('EPS+1', (EPS + 1) != 1, True)

        self.test('EPS1',    EPS1 < 1, True)
        self.test('EPS1-1', (EPS1 - 1) != 1, True)

        self.test('.tillC', _AtoZnoIO_.tillC, 'ABC')
        self.test('.fromX', _AtoZnoIO_.fromX, 'XYZ')

        self.test('.fromH.tillJ', _AtoZnoIO_.fromH.tillJ, 'HJ')
        self.test('.fromN.tillP', _AtoZnoIO_.fromN.tillP, 'NP')


if __name__ == '__main__':

    t = Tests(__file__, __version__, interns)
    t.testInterns()
    t.results()
    t.exit()
