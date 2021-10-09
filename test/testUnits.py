
# -*- coding: utf-8 -*-

# Test units module.

__all__ = ('Tests',)
__version__ = '21.09.30'

from base import TestsBase

from pygeodesy import Band, Bearing, Bearing_, Bool, \
                      Epoch, Epsg, FIx, Garef, Geohash, Georef, \
                      Int, Int_, Number_, Precision_, Lam_, Phi_, \
                      Str, Zone, Float, units
_NamedUnit = units._NamedUnit


class Tests(TestsBase):

    def testUnit(self, U, arg, known=False):
        self.subtitle(units, 'ing %s%s ' % (U.__name__, (arg,)))

        r = repr(arg)
        s = str(arg)

        n = U.__name__.lower()
        u = U(arg, name=n)
        u.units = n
        R = '%s (%s)' % (u.name, r)

        self.test('.classname', u.classname, U.__name__)
        self.test('isinstance', isinstance(u, U), True)
        self.test('.name', u.name, n)
        self.test('.named', u.named, n)
        self.test('.named2', u.named2, u.classname + ' %r' % (n,))
        self.test('.str', str(u), s)
        self.test('.toStr', u.toStr(), s)
        self.test('.std_repr', u.std_repr, u.std_repr)
        self.test('.repr', repr(u), r if u.std_repr else R)
        self.test('.toRepr', u.toRepr(), R)
        self.test('.units', u.units, n)

        try:
            self.test('error', U('X'), ValueError, known=known)
        except ValueError as x:
            self.test(n, str(x), str(x))
        try:
            self.test('error', U('X', name='U'), ValueError, known=known)
        except ValueError as x:
            self.test(n, str(x), str(x))
        try:
            self.test('Error', U('X', Error=TypeError), TypeError, known=known)
        except TypeError as x:
            self.test(n, x.__class__.__name__, TypeError.__name__)

        u.rename('Test')
        R = '%s (%r)' % (u.name, arg)

        self.test('.named', u.named, 'Test')
        self.test('.named2', u.named2, u.classname + ' %r' % ('Test',))
        self.test('.str', str(u), s)
        self.test('.toStr', u.toStr(), s)
        self.test('.repr', repr(u), r if u.std_repr else R)
        self.test('.toRepr', u.toRepr(), R)
        self.test('.units', u.units, n)

        for a in ('name', '_name'):  # coverage
            self.test('.'+a, getattr(u, a), 'Test')
        try:  # coverage
            self.test('.str', str(u), u)
        except AssertionError as x:
            self.test('.str', x, x)  # PYCHOK test
        try:  # coverage
            self.test('.repr', repr(u), r if u.std_repr else R)
        except AssertionError as x:
            self.test('.repr', x, x)

        delattr(u, '_name')  # coverage
        self.test('delattr', repr(u.name), "''")

    def testUnits(self):
        for U in self.pygeodesy_classes(_NamedUnit):
            if U not in (Band, Bool, Bearing_,
                         Epoch, Epsg, FIx,
                         Garef, Geohash, Georef,
                         Int, Int_, Number_,
                         Precision_, Str,
                         Lam_, Phi_, Zone,
                        _NamedUnit):
                self.testUnit(U, 1.0)  # sample

        for U in (Band, Str):
            self.testUnit(U, 'U', known=True)

        for U in (Bool,):
            self.testUnit(U, True, known=True)

        for U in (Int, Int_, Number_, Precision_, Zone):
            self.testUnit(U, 2)

        for U in (Epoch,):
            self.testUnit(U, 1901, known=True)

        self.subtitle(units)  # courtesy JaapZee at Gmail
        self.test(Bearing.__name__,  Bearing(361), 1.0)
        self.test(Bearing_.__name__, Bearing_(361), 0.01745, fmt='%.5f')

        self.test(Lam_.__name__, Lam_(361, clip=0), 6.3, fmt='%.2f')
        self.test(Phi_.__name__, Phi_(361, clip=0), 6.3, fmt='%.2f')

        self.test(FIx.__name__, FIx(1),   Int(1), known=True)
        self.test(FIx.__name__, FIx(1.5), Float(1.5),)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testUnits()
    t.results()
    t.exit()
