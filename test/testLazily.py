
# -*- coding: utf-8 -*-

# Test the lazy import module L{lazily}.

__all__ = ('Tests',)
__version__ = '23.11.30'

from bases import TestsBase, _env_c2, isPython37, type2str
import pygeodesy


class Tests(TestsBase):

    def testLazily(self):

        _all_  = pygeodesy.__all__
        lazily = pygeodesy.lazily

        for a in sorted(_all_, key=str.lower):
            t = type2str(pygeodesy, a).replace('()', '').strip()
            self.test(a, t, t)

        z = pygeodesy.isLazy
        self.test('isLazy', z, z, nl=1)
        if not z:
            for a, m in lazily._all_missing2(_all_):
                t = all(_.startswith('rhumb') for _ in m) if m else False
                m = ', '.join(m) if m else None
                self.test('missing in %s' % (a,), m, None, known=t)

        # simple lazy_import enable tests
        env_cmd, cmd = _env_c2("-c 'import sys; import pygeodesy; "
                                   "sys.exit(0 if pygeodesy.isLazy == %s else 1)'")
        self.test('cmd', cmd, cmd, nl=1)
        if env_cmd:
            from os import system
            for z in range(5):
                e = 'PYGEODESY_LAZY_IMPORT=%s' % (z,)
                c =  env_cmd % (e, z if isPython37 else None)
                x =  system(c) >> 8  # exit status in high byte
                self.test(e, x, 0)
        else:
            self.skip('no env_cmd')

        M = lazily._ALL_MODS
        self.test('items', M.name, M.name, nl=1)
        for n, m in sorted(M.items()):  # coverage
            self.test(n, m, m)

#       p = pygeodesy.printf
#       t = p('to %(std)s, flushed', nl=1, std='stdout', file=sys.stdout, flush=True)
#       self.test('%s to stdout, flushed' % (p.__name__,), t, t)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testLazily()
    t.results()
    t.exit()
