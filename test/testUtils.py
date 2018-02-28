
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '18.02.27'

from base import TestsBase

from pygeodesy import R_M, antipode, fStr, heightOf, horizon, \
                         isantipode, unroll180


class Tests(TestsBase):

    def testUtils(self):

        self.test('antipode1', antipode( 89,  179), (-89, -1))
        self.test('antipode2', antipode(-89, -179),  (89,  1))

        self.test('isantipode1', isantipode( 89,  179, -89, -1), True)
        self.test('isantipode2', isantipode(-89, -179,  89,  1), True)
        self.test('isantipode3', isantipode(-89, -179, -89, -1), False)

        self.test('heightof0',   heightOf(0,   R_M), 2638958.23912, fmt='%.5f')
        self.test('heightof45',  heightOf(45,  R_M), 5401080.43931, fmt='%.5f')
        self.test('heightof90',  heightOf(90,  R_M), R_M)
        self.test('heightof135', heightOf(135, R_M), 5401080.43931, fmt='%.5f')

        self.test('horizon0',     horizon(0), 0.0)
        self.test('horizon10Km',  horizon(10000), '357099.672', fmt='%.3f')
        self.test('horizon30Kft', horizon(10000, refraction=True), '392310.704', fmt='%.3f')
        self.test('horizon10Kft', horizon( 3000, refraction=True), '214877.422', fmt='%.3f')

        self.test('unroll180', fStr(unroll180(-90, 110, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 110, wrap=False)), '200.0, 110.0')

        self.test('unroll180', fStr(unroll180(-90, 830, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 830, wrap=False)), '920.0, 830.0')

        self.test('unroll180', fStr(unroll180(-110, 90, wrap=True)), '-160.0, -270.0')
        self.test('unroll180', fStr(unroll180(-110, 90, wrap=False)), '200.0, 90.0')

        self.test('unroll180', fStr(unroll180(-830, 90, wrap=True)), '-160.0, -990.0')
        self.test('unroll180', fStr(unroll180(-830, 90, wrap=False)), '920.0, 90.0')


if __name__ == '__main__':

    from pygeodesy import utils  # private

    t = Tests(__file__, __version__, utils)
    t.testUtils()
    t.results(nl=0)
    t.exit()
