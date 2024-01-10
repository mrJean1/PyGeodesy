
# -*- coding: utf-8 -*-

# Test _isfrozen for PyInstaller, etc. errors.

__all__ = ('Tests',)
__version__ = '23.12.30'

from bases import TestsBase, _env_c2


class Tests(TestsBase):

    def testFrozen(self):
        # simple _isfrozen test, mimicking PyInstaller
        env_cmd, cmd = _env_c2("-c 'import sys; sys.frozen = True; import pygeodesy; "
                                   "sys.exit(0 if pygeodesy._isfrozen else 1)'")
        self.test('cmd', cmd, cmd)
        if env_cmd:
            from os import system
            for z in range(3):
                e = 'PYGEODESY_LAZY_IMPORT=%s' % (z,)
                c =  env_cmd % (e,)
                x =  system(c) >> 8  # exit status in high byte
                self.test(e, x, 0)
        else:
            self.skip('no env_cmd')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testFrozen()
    t.results()
    t.exit()
