
# -*- coding: utf-8 -*-

u'''Test properties.
'''

__all__ = ('Tests',)
__version__ = '21.01.08'

from base import TestsBase

from pygeodesy import Property, Property_RO, property_RO, props


class Tests(TestsBase):

    def testProps(self, Base, *args):  # MCCABE 26

        self.subtitle(props, Base)

        class T(Base):

            _p = 0
            _q = 1

            @Property_RO
            def P(self):
                return self.p  # PYCHOK property

            @property_RO
            def p(self):
                self._p += 1  # PYCHOK _p?
                return self._p  # PYCHOK _p?

            @Property
            def q(self):
                return self._q  # PYCHOK _q?

            @q.setter  # PYCHOK .setter
            def q(self, q):
                self._q = q

        t = T(*args)
        self.test('P1', t.P, 1)
        self.test('p1', t.p, 2)

        self.test('P2', t.P, 1)
        self.test('p2', t.p, 3)

        self.test('q1', t.q, 1)
        t.q = 2
        self.test('q2', t.q, 2)
        self.test('q3', 'q' in t.__dict__, True)

        T.P._update(t)
        self.test('P3', t.P, 4)
        self.test('p3', t.p, 5)

        T.q._update(t)
        self.test('q4', 'q' in t.__dict__, False)
        self.test('q5', t.q, 2)

        class E(Base):

            @Property_RO
            def X(self):
                return None

            try:
                @X.setter  # PYCHOK .setter
                def X(self, x):
                    self._x = x

                self.test('X1', 'X.setter', AttributeError)
            except AttributeError as x:
                self.test('X1', str(x), 'immutable Property_RO: X.setter X')  # PYCHOK x?

            try:
                @X.deleter  # PYCHOK .deleter
                def X(self):
                    pass

                self.test('X2', 'X.deleter', AttributeError)
            except AttributeError as x:
                self.test('X2', str(x), 'invalid Property_RO: X.deleter X')

            @property_RO
            def y(self):
                return None

            try:
                @y.setter  # PYCHOK .setter
                def y(self, y):
                    self._y = y

                self.test('y1', 'y.setter', AttributeError)
            except AttributeError as x:
                self.test('y1', str(x), 'immutable property_RO: y.setter y')

            try:
                @y.deleter  # PYCHOK .deleter
                def y(self):
                    pass

                self.test('y2', 'y.deleter', AttributeError)
            except AttributeError as x:
                self.test('y2', str(x), 'invalid property_RO: y.deleter y')

            @Property
            def Z(self):
                return None

            try:
                @Z.deleter  # PYCHOK .deleter
                def Z(self):
                    pass

                self.test('Z1', 'Z.deleter', AttributeError)
            except AttributeError as x:
                self.test('Z1', str(x), 'invalid Property: Z.deleter Z')

            try:
                @Z.getter  # PYCHOK .getter
                def Z(self):
                    pass

                self.test('Z2', 'Z.getter', AttributeError)
            except AttributeError as x:
                self.test('Z2', str(x), 'invalid Property: Z.getter Z')

        del E


if __name__ == '__main__':

    from pygeodesy import Ellipsoid, R_M
    from pygeodesy.named import _NamedBase

    t = Tests(__file__, __version__)
    t.testProps(_NamedBase)
    t.testProps(Ellipsoid, R_M, R_M)
    t.results()
    t.exit()
