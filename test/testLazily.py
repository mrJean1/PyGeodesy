
# -*- coding: utf-8 -*-

# Test the lazy import module lazily.

__all__ = ('Tests',)
__version__ = '20.01.07'

from base import TestsBase, ismacOS, isNix, isPython37, isWindows, \
                 PythonX, type2str
import pygeodesy

lazily = pygeodesy.lazily
_all_  = pygeodesy.__all__

import os

_cmd = PythonX + " -c 'import pygeodesy, sys; " \
                      "sys.exit(0 if pygeodesy.isLazy == %s else 1)'"
if ismacOS or isNix:
    _env_cmd = 'env %s ' + _cmd + ' >>/dev/null'
elif isWindows:  # XXX UNTESTED
    _env_cmd = 'set %s;' + _cmd
else:
    _env_cmd = None

_HOME = os.environ.get('HOME', '')
if _HOME and _cmd.startswith(_HOME):
    _cmd = '~' + _cmd[len(_HOME):]
del _HOME


class Tests(TestsBase):

    def testLazily(self):

        for a in sorted(_all_, key=str.lower):
            t = type2str(pygeodesy, a).replace('()', '').strip()
            self.test(a, t, t)

        z = pygeodesy.isLazy
        self.test('isLazy', z, z)
        if not z:
            for a, m in lazily._all_missing2(_all_):
                self.test('missing in %s' % (a,), m or None, None)

        # simple lazy_import enable tests
        self.test('cmd', _cmd, _cmd)
        if _env_cmd:
            for z in range(5):
                e = 'PYGEODESY_LAZY_IMPORT=%s' % (z,)
                c = _env_cmd % (e, z if isPython37 else None)
                x = os.system(c) >> 8  # exit status in high byte
                self.test(e, x, 0)
        else:
            self.skip('no _env_cmd')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testLazily()
    t.results()
    t.exit()
