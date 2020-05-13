
# -*- coding: utf-8 -*-

# Test formatting functions.

__all__ = ('Tests',)
__version__ = '20.05.03'

from base import TestsBase

from pygeodesy import INF, NEG0, NAN, \
                      anstr, fstr, fstrzs, instr, LatLon_, unstr


class Tests(TestsBase):

    def testStreprs(self):

        self.test('anstr', anstr('a-b?_'), 'a-b__')

        self.test('fstr', fstr(0.123, prec=-6), '0.123000')
        self.test('fstr', fstr(0.123, prec=+6), '0.123')
        self.test('fstr', fstr((0.123, 456.789), prec=+6), '0.123, 456.789')
        self.test('fstr', fstr(0.123, prec=-5, fmt='%.*e'), '1.23000e-01')
        self.test('fstr', fstr(0.123, prec=+5, fmt='%.*e'), '1.23e-01')
        try:  # coverage
            self.test('fstr', fstr(1, fmt='X'), ValueError.__name__)
        except ValueError as x:
            self.test('fstr', str(x), "fmt ('X'): not '[%.*]F|f|E|e|G|g'")

        for f, x in ((1,      '1.0'),
                     (1.0,    '1.0'),
                     (-1,    '-1.0'),
                   # (1e300,  '10000<290>40160.'),
                     (1e1234, 'INF'),  # == INF
                     (INF,    'INF'),
                     (NAN,    'NAN'),
                     (NEG0,  '-0.0'),
                     (0.0,    '0.0')):
            self.test('fstr(%F)' % (f,), fstr(f), x)

        for f, x in (('0.0',      '0.0'),
                     ('0.00',     '0.0'),
                     ('0.000',    '0.0'),
                     ('00.0',    '00.0'),
                     ('000.00', '000.0'),
                     ('0.000',    '0.0'),
                     ('0.010',    '0.01'),
                     ('0.0200',   '0.02'),
                     ('0.0e+01',      '0.0e+01'),
                     ('0.00e+02',     '0.0e+02'),
                     ('0.000e+03',    '0.0e+03'),
                     ('00.0e+00',    '00.0e+00'),
                     ('000.00e+01', '000.0e+01'),
                     ('0.000e+02',    '0.0e+02'),
                     ('0.010e+03',    '0.01e+03'),
                     ('0.0200e+00',   '0.02e+00')):
            self.test('fstrzs(%s)' % (f,), fstrzs(f), x)

        for f, x in (('0',    '0.0'),
                     ('0.0',  '0.0'),
                     ('0.',   '0.'),
                     ('1e10', '1.0e10'),
                     ('2E+2', '2.0E+2'),
                     ('3.E3', '3.E3')):
            self.test('fstrzs(%s, ap1z=True)' % (f,), fstrzs(f, ap1z=True), x)

        ll = LatLon_(45, 90, height=1.2)
        self.test('instr', ll.toStr2(), 'LatLon_(45.0°N, 090.0°E, +1.20)')
        self.test('instr', instr(ll, 45, 90, h=1.2), 'LatLon_(45, 90, h=1.2)')

        self.test('unstr', unstr('f', 1.1, 2.2), 'f(1.1, 2.2)')
        self.test('unstr', unstr('f', y=2.2, x=1.1), 'f(x=1.1, y=2.2)')


if __name__ == '__main__':

    from pygeodesy import streprs

    t = Tests(__file__, __version__, streprs)
    t.testStreprs()
    t.results()
    t.exit()
