
# -*- coding: utf-8 -*-

u'''Test L{interns} module.
'''

__all__ = ('Tests',)
__version__ = '22.04.21'

from base import TestsBase

from pygeodesy import clips, interns, isinf, isnan, EPS, EPS0, EPS02, \
                      EPS1, EPS2, EPS_2, EPS4, INF, INT0, NAN, NEG0, NINF, NN

from os import getcwd

_0to9_        =  interns._0to9_
_AtoZnoIO_    =  interns._AtoZnoIO_
_cwd          =  getcwd()
_DUNDER_      =  interns._DUNDER_
_DOT_         =  interns._DOT_
_EQUALSPACED_ =  interns._EQUALSPACED_
_functions    = (interns._90_EPS_2,
                 interns._float,
                 interns._float0,
                 interns._floatuple,
                 interns._platform2,
                 interns._version2)
_UNDER_       =  interns._UNDER_
_exceptions   = (_0to9_, _AtoZnoIO_,
                 interns._doesn_t_exist_,
                 interns._exceed_PI_radians_,
                 interns._iadd_,
                 interns._n_a_,
                 interns._NL_hash_,
                 interns._NL_var_,
                 interns._not_finite_,
                 interns._not_scalar_,
                 interns._PyPy__,
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
                self.test(n, isinstance(i, (float, str)) or i in _functions, True, known=k)
                # check the naming conventions
                if isinstance(i, str) and not k:
                    a = a.strip(_UNDER_)
                    self.test(n, i.lower(), a.lower(), known=a.isupper() or i in _exceptions)
                elif isinstance(i, float):
                    a = a.strip(_UNDER_).replace(_UNDER_, _DOT_)
                    self.test(n, i, a, known=n.startswith(_UNDER_))

        self.test('EPS',    EPS > 0, True)
        self.test('EPS+1', (EPS + 1) != 1, True)

        self.test('EPS0',  0 < EPS0  < EPS,  True)
        self.test('EPS02', 0 < EPS02 < EPS0, True)

        self.test('EPS_2', 0 < EPS_2 < EPS, True)
        self.test('EPS_2',     EPS_2,  EPS / 2)

        self.test('EPS1',    EPS1 < 1, True)
        self.test('EPS1-1', (EPS1 - 1) != 1, True)

        self.test('EPS2', EPS  < EPS2 < EPS4, True)
        self.test('EPS4', EPS2 < EPS4, True)

        self.test('INF',  isinf(INF),  True)
        self.test('INF',  INF == NINF, False)
        self.test('NINF', isinf(NINF), True)
        self.test('NINF',       NINF, -INF)

        self.test('INT0', INT0, '0')

        self.test('NAN',  isnan(NAN),  True)
        self.test('NAN',  NAN == INF,  False)
        self.test('NAN',  NAN == NINF, False)

        self.test('NEG0',       NEG0,  '-0.0')
        self.test('NEG0',  0 == NEG0,   True)
        self.test('NEG0',  bool(NEG0),  False)
        self.test('NEG0',   abs(NEG0), '0.0')

        self.test('.tillC', _AtoZnoIO_.tillC, 'ABC')
        self.test('.fromX', _AtoZnoIO_.fromX, 'XYZ')

        self.test('.fromH.tillJ', _AtoZnoIO_.fromH.tillJ, 'HJ')
        self.test('.fromN.tillP', _AtoZnoIO_.fromN.tillP, 'NP')

        _90_EPS_2 = interns._90_EPS_2
        self.test(_90_EPS_2.__name__, _90_EPS_2(90) < 90, True)
        self.test(_90_EPS_2.__name__, _90_EPS_2(90) > 89.999999, True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, interns)
    t.testInterns()
    t.results()
    t.exit()
