
# -*- coding: utf-8 -*-

# Test units module.

__all__ = ('Tests',)
__version__ = '20.05.10'

from base import TestsBase

from pygeodesy import (Bearing, Bearing_,
                       Degrees, Distance, Easting,
                       Feet, Float, Height, Int,
                       Lam, Lam_, Lat, Lon, Meter, Northing,  # Lam_
                       Phi, Phi_, Radians, Radius, Radius_,  # Phi_
                       Scalar, Scalar_, Number_, Precision_,
                       Str, units)


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
        self.test('.repr', repr(u), r if u.std_repr else R)
        self.test('.toRepr', u.toRepr(), R)
        self.test('.units', u.units, n)

        try:
            self.test('error', U('X'), ValueError, known=known)
        except ValueError as x:
            self.test(n, str(x), str(x))  # PYCHOK test
        try:
            self.test('error', U('X', name='U'), ValueError, known=known)
        except ValueError as x:
            self.test(n, str(x), str(x))  # PYCHOK test
        try:
            self.test('Error', U('X', Error=TypeError), TypeError, known=known)
        except TypeError as x:
            self.test(n, x.__class__.__name__, TypeError.__name__)  # PYCHOK test

        u.name = 'Test'
        R = '%s (%r)' % (u.name, arg)

        self.test('.named', u.named, 'Test')
        self.test('.named2', u.named2, u.classname + ' %r' % ('Test',))
        self.test('.str', str(u), s)
        self.test('.toStr', u.toStr(), s)
        self.test('.repr', repr(u), r if u.std_repr else R)
        self.test('.toRepr', u.toRepr(), R)
        self.test('.units', u.units, n)

        for a in ('name', '_name'):  # coverage
            self.test('.'+a, getattr(u, a), 'Test')  # PYCHOK test
        try:  # coverage
            self.test('.str', str(u), u)
        except AssertionError as x:
            self.test('.str', x, x)  # PYCHOK test
        try:  # coverage
            self.test('.repr', repr(u), r if u.std_repr else R)  # PYCHOK test
        except AssertionError as x:
            self.test('.repr', x, x)

        delattr(u, '_name')  # coverage
        self.test('delattr', repr(u.name), "''")

    def testUnits(self):
        for U in (Float, Bearing,
                         Degrees, Distance, Easting, Feet, Height,  # PYCHOK indent
                         Lam, Lat, Lon, Meter, Northing,  # PYCHOK indent
                         Phi, Radians, Radius, Radius_,  # PYCHOK indent
                         Scalar, Scalar_):  # PYCHOK indent
            self.testUnit(U, 1.0)

        for U in (Int, Number_, Precision_):
            self.testUnit(U, 2)

        for U in (Str,):
            self.testUnit(U, 'U', known=True)

        self.subtitle(units)  # courtesy JaapZee at Gmail
        self.test(Bearing.__name__,  Bearing(361), 1.0)
        self.test(Bearing_.__name__, Bearing_(361), 0.01745, fmt='%.5f')
        self.test(Lam_.__name__,     Lam_(361, clip=0), 6.3, fmt='%.1f')
        self.test(Phi_.__name__,     Phi_(361, clip=0), 6.3, fmt='%.1f')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testUnits()
    t.results()
    t.exit()
