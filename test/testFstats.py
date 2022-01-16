
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.01.14'

from base import TestsBase

from pygeodesy import EPS0, Fcook, Flinear, fstats, Fsum, Fwelford


class Tests(TestsBase):

    def testFstats(self):

        c = Fcook((2, 4, 4, 4, 5, 5, 7, Fsum(9)), name=Fcook.__name__)
        self.test(c.name, len(c), 8)
        k = c.fkurtosis()
        self.test(c.name, k, '-0.218750', prec=6)
        e = c.fmedian()
        self.test(c.name, e, '4.562500', prec=6)
        m = c.fmean()
        self.test(c.name, m, '5.0', known=abs(m - 5) < EPS0)
        s = c.fskewness()
        self.test(c.name, s, '0.656250', prec=6)
        d = c.fstdev(sample=False)
        self.test(c.name, d, '2.0', prec=1)
        v = c.fvariance()
        self.test(c.name, v, d**2, prec=1)
        j = c.fjb()
        self.test(c.name, j, '1.039635', prec=6, nt=1)

        t = c.fcopy(name=c.fcopy.__name__)
        self.test(t.name, t, c, known=True)
        self.test(t.name, t.fmean_(),  m)
        self.test(t.name, t.fstdev_(), d, prec=1)
        self.test(c.name, c.fadd(()), 8, nt=1)

        c = t + c
        c.rename('Doubled')
        self.test(c.name, len(c),  16)
        self.test(c.name, c.fkurtosis(), k, prec=6)
        self.test(c.name, c.fmedian(), e)
        self.test(c.name, c.fmean(), m)
        self.test(c.name, c.fskewness(), s, prec=6)
        self.test(c.name, c.fstdev(), d, prec=1)
        self.test(c.name, c.fvariance(), v, prec=1)
        self.test(c.name, c.fjb(), j, prec=6, known=True, nt=1)

        t = Fcook(name='Empty')
        t += c
        self.test(t.name, t, c, known=True)
        t += 0
        t += Fsum()
        self.test(t.name, len(t), 18)
        try:
            t += None
            self.test(t.name, t, TypeError.__name__)
        except TypeError as x:
            self.test(t.name, str(x), str(x))

        # <https://www.Real-Statistics.com/descriptive-statistics/symmetry-skewness-kurtosis/>
        c = Fcook((2, 5, -1, 3, 4, 5, 0, Fsum(2)), name='Excel')
        self.test(c.name, len(c), 8, nl=1)
        k = c.fkurtosis_()
        self.test(c.name, k, '-1.114187', prec=6)
        t = c.fkurtosis_(sample=True)
        self.test(c.name, t, '-0.939792', prec=6)
        e = c.fmedian_()
        self.test(c.name, e, '2.735294', prec=6)
        m = c.fmean_()
        self.test(c.name, m, '2.50', prec=2, known=abs(m - 5) < EPS0)
        s = c.fskewness_()
        self.test(c.name, s, '-0.342403', prec=6)
        t = c.fskewness_(sample=True)
        self.test(c.name, t, '-0.427052', prec=6)
        d = c.fstdev_(sample=False)
        self.test(c.name, d, '2.061553', prec=6)
        v = c.fvariance_()
        self.test(c.name, v, d**2, prec=1)
        j = c.fjb_()
        self.test(c.name, j, '0.470372', prec=6, nt=1)

        w = Fwelford((2, 4, 4, 4, 5, 5, 7, Fsum(9)), name=Fwelford.__name__)
        self.test(w.name, len(w), 8)
        m = w.fmean()
        self.test(w.name, m, '5.0', known=abs(m - 5) < EPS0)
        d = w.fstdev(sample=False)
        self.test(w.name, d, '2.0', prec=1, known=abs(m - 2) < EPS0)
        v = w.fvariance(sample=False)
        self.test(w.name, v, d**2, prec=1, nt=1)

        t = w.fcopy(name=w.fcopy.__name__)
        self.test(t.name, t, w, known=True)
        self.test(t.name, t.fmean_(),  m)
        self.test(t.name, t.fstdev_(), d, prec=1)
        self.test(t.name, t.fvariance_(), d**2, prec=1)
        self.test(t.name, t.fadd_(), 8)

        t = Fwelford(name='Empty') + c.toFwelford()
        t += t
        t += 0
        t += Fsum()
        self.test(t.name, len(t), 18)
        try:
            t += None
            self.test(t.name, t, TypeError.__name__)
        except TypeError as x:
            self.test(t.name, str(x), str(x))

        # <https://www.DisplayR.com/what-is-linear-regression/>
        r = Flinear(( 23,  26,  30,   34,   43,   48,   52,   57,   58),
                    (651, 762, 856, 1063, 1190, 1298, 1421, 1440, 1518), name=Flinear.__name__)
        self.test(r.name, len(r), 9, nl=1)
        self.test(r.name, r.fcorrelation(), '0.988288', prec=6)
        self.test(r.name, r.fintercept(), '167.682949', prec=6)
        self.test(r.name, r.fslope(),      '23.422786', prec=6)

        r = r + (1, 2)
        self.test(r.name, len(r), 10)
        r += r.fcopy()
        r += Fsum(), Fsum()
        self.test(r.name, len(r), 21)
        try:
            r += None
            self.test(r.name, r, TypeError.__name__)
        except TypeError as x:
            self.test(r.name, str(x), str(x))
        try:
            r += 1, 2, 3
            self.test(r.name, r, ValueError.__name__)
        except ValueError as x:
            self.test(r.name, str(x), str(x))


if __name__ == '__main__':

    t = Tests(__file__, __version__, fstats)
    t.testFstats()
    t.results()
    t.exit()
