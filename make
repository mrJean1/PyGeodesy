#!/usr/bin/env python3.13
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests, rebuild the docs, etc.

from collections import defaultdict
from contextlib import contextmanager
from glob import iglob
import os
from subprocess import PIPE, STDOUT, Popen
from time import sleep, time
import sys

__version__ = '24.10.20'
_PyGeodesy_ = 'PyGeodesy'

_OO = '' if __debug__ else ' -OO'
_PD = os.getcwd().split(_PyGeodesy_)[0]
_SIsecs = 'fs', 'ps', 'ns', 'us', 'ms', 'sec'  # reversed
_time0  = time()


def _all_locals(*files):  # PYCHOK no cover
    # check __all__ vs locals().
    from os.path import basename, splitext

    import pygeodesy  # Python 2- only
    if not pygeodesy.lazily._FOR_DOCS:
        sys.exit('usage: env PYGEODESY_FOR_DOCS=1 make ...')

    for f in (files or _walker(1)):  # sorted?
        n = splitext(basename(f))[0]
        m = getattr(pygeodesy, n)
        d = dict((a, getattr(m, a)) for a in dir(m))
        for a in m.__all__:
            if a not in d and not a.startswith('_'):
                _printf('%s.%s %r not in %s.%s', n,'__all__', a, n,'locals')
        for a, o in d.items():
            if a not in m.__all__ and getattr(o, '__module__', '') == n:
                _printf('%s.%s %r not in %s.%s', n,'locals', a, n,'__all__')


def _cmd(cmd, *args):
    if args:
        cmd = cmd % args
    _gc(cmd)
    return os.system(cmd)


def _coverage_percentage(html, readme):
    # check coverage percentage in README.rst
    p = ''
    t = '<span class="pc_cov">'
    h = _read(html)
    i =  h.find(t)
    if i > 0:
        i += len(t)
        j = h.find('%</span>', i)
        if j > i:
            p = h[i:j].strip() + '%'
    if not p:
        sys.exit('not in %s: %r' % (html, t))

    r = _read(readme)
    t = 'https://Img.Shields.io/badge/coverage-%s25-brightgreen' % (p,)
    if t not in r:
        sys.exit('not in %s: %r' % (readme, t))

    _cmd('echo coverage: %s %r in %s OK', html, p, readme)


def _coverage_pygeodesy(*cb_, **z_p_):  # MCCABE 15
    # remove _cb_...ext tails
    cb = {}
    for c in cb_:
        i = c.find('_cb_')
        j = c.rfind('.')
        if 0 < i < j < len(c):
            b = c[:i] + c[j:]
            if b != c:
                os.system('mv testcoverage/' + c +
                            ' testcoverage/' + b)
                cb[c] = b
    # replace 'z_4569faf12939b165_' with 'pygeodesy_'
    os.system('rm -f testcoverage/class_index.html'
                   ' testcoverage/function_index.html')
#   for h in iglob('testcoverage/*index.html'):
#       t = r = _read(h)
#       for c, b in cb.items():
#           t = t.replace(c, b)
#       _write(h, t, r)

    for h in iglob('testcoverage/*.html'):
        for z_, p_ in z_p_.items():
            t = h.replace(z_, p_)
            if t != h:
                os.system('mv ' + h + ' ' + t)
                h = t
                r = t = _read(h)
                t = t.replace(z_, p_)
                for c, b in cb.items():
                    t = t.replace(c, b)
                _write(h, t, r)
                break

    h = 'testcoverage/index.html'
    r = t = _read(h)
    i = t.find('<h2>')
    if i > 0:  # buttons Classes, etc.
        j = t.find('</h2>', i) + 5
        if j > i:
            t = t[:i] + t[j:]
    for z_, p_ in z_p_.items():
        t = t.replace(z_, p_)
    for c, b in cb.items():
        t = t.replace(c, b)
    _write(h, t, r)

    # take out <aside ...> to </form> tags
    i = t.find('<aside ')
    if i > 0:
        j = t.find('</form>', i) + 7
        if j > i:
            t = t[:i] + t[j:]
    return t


def _coverage_pygeodesy_cleaner():
    x = _coverage_pygeodesy('coverage_html_cb_6fb7b396.js', 'favicon_32_cb_58284776.png',
                            'keybd_closed_cb_ce680311.png', 'style_cb_8e611ae1.css',
                            z_4569faf12939b165_='pygeodesy_',
                            z_5087a5e01bbcffc5_='pygeodesy_auxilats_',
                            z_81bb61ef01d6cb97_='pygeodesy_deprecated_',
                            z_4d2bca0ead806276_='pygeodesy_geodesicx_',
                            z_5cf368818b86873a_='pygeodesy_rhumb_')  # coverage 7.2.2+ only
    return x


def _dist(dirname):
    # check that all sub=packackes are included, everywhere needed
    for n, sps, _ in os.walk(dirname):
        break
    if n != dirname:
        sys.exit('package %r bad: %r' % (dirname, n))
    try:
        sps.remove('__pycache__')
    except ValueError:
        pass

    from pygeodesy.interns import _SUB_PACKAGES
    if tuple(sorted(sps)) != tuple(sorted(_SUB_PACKAGES)):
        sys.exit('%s: lazily %r vs %s %r' % (_dist.__name__, _SUB_PACKAGES, n, sps))

    _init_py = '%s/__init__.py' % (n,)
    for f, x in (('setup.py',        ", '%s.%s'"),
                 ('MANIFEST.in',     'graft %s/%s'),
                 ('testcoverage.rc', '    %s/%s/'),
                 (_init_py,          '    import %s.%s'),
                 (_init_py,          '   from %s.%s ')):  # auxilats
        t = _read(f)
        for p in sps:
            p = x % (n, p)
            if p not in t:
                sys.exit('%s: %r missing in %s' % (_dist.__name__, p, f))


def _eLs(py, pubs):
    # all public L{...} exports
    for n, s, t in _readlines3(py):
        j = 0
        while True:
            i = s.find('L{', j)
            if i < 0:
                break
            i += 2
            j  = s.find('}', i)
            if j > i:
                L = s[i:j]
                p, _, g = L.partition('.')
                if not g:
                    g = p
                elif p != 'pygeodesy':
                    g = L
                if g in pubs or ('def %s(' % (L,)) in t or \
                                ('class %s(' % (L,)) in t:
                    continue
                yield n, L, s
            else:
                j = i


def _epytext(py, d, *Ls):
    # check for balanced {...} brackets
    t = _read(py)
    j = 0
    while True:
        i = t.find('{', j) - 1
        if i < 0:
            break
        if t[i] in Ls:  # 'BCILMU':
            j = t.find('}', i + 2)
            if j > i + 1:
                if '{' in t[i + 2:j]:
                    j = t.find('}', j + 1, j + 100)
            if j > i:
                m = ' '.join(t[i:j + 1].split())
                if m[-1] != '}':
                    m += ' ***'
                n = ':' + str(t[:i].count('\n') + 1)
                if m in d:
                    d[m].append(py + n)
                else:
                    d[m] = [py + n]
            else:
                j = i + 2
        else:
            j = i + 2
    return d


def _files(files):
    d = {}
    for f in files:
        d[f] = d.get(f, 0) + 1
    return '  '.join(_filesort(d))


def _filesort(d):
    for f, c in sorted(d.items()):
        if c > 1:
            f += ' %s' % (c,)
        yield f


def _froms():
    # check from <module> import ...
    d = defaultdict(tuple)
    for p in _walker(1):
        for n, s, _ in _readlines3(p):
            s = s.strip().split(None, 4)
            if (len(s) >= 2 and s[0] == 'import' or (
                len(s) >= 3 and s[0] == 'from'
                            and s[2] == 'import')):
                d[s[1]] += (p + ':' + str(n),)

    from importlib import import_module
    for m, t in sorted(d.items()):
        try:
            i = import_module(m)
            del i
        except ImportError as x:
            print('%s: from %s import ... in %s' % (x, m, ' '.join(t)))


def _gc(cmd):
    if len(cmd) > 104:
        cmd = cmd[:50] + '...' + cmd[-50:]
    _printf('\ncleaning for%s: %s (%s) ...', _OO, cmd, _secs2str(time() - _time0))
    c = '*.o *.pyc *.pyo *.*~ .DS_Store ._* .*ignore'.split()
    for g in (' ', ' */', ' */*/'):
        os.system('rm -f ' + g.join(c))
        os.system('rm -rf %s__pycache__' % (g,))


def _html_cut(html, tag, m, end, n):
    # cut off some tags
    j = 0
    while True:
        i = html.find(tag, j)
        if i < 0:
            break
        i += m
        j  = html.find(end, i, i+128)
        if j < 0:
            j = i
        else:
            j += n
            html = html[:i] + html[j:]
    return html


_html_from_tos = (
    ('<i>unreachable</i>.', ''),
    (_PD, '..../'),
    ('<i>attributes-LatLon-html</i>', '<a href="attributes-LatLon.html"><i>LatLon</i></a>'),
    ('>\\xc2\\xb0<',         '>&deg;<'),
    ('>\\xe2\\x80\\xb3<',    '>&Prime;<'),  # &#8243;
    ('<dt>Parameters:</dt>', '<dt>Arguments:</dt>'),
    ('<dt>Get Method:</dt>', '<dt>Get method:</dt>'),
    ('<dt>Set Method:</dt>', '<dt>Set method:</dt>'),
    ('</span>(<span class="sig-arg">', '</span>&nbsp;(<span class="sig-arg">'))


def _html(html, *from_tos):
    # tweak html file
    t =  h = _read(html)
    t = _html_cut(t, '</a>.  &lt;function ', 5, '&gt;</p>', 4)
    for r in from_tos:
        t = t.replace(*r)
    if t != h:
        _printf('editing: %s', html)
        _write(html, t)


def _long_description():
    t = ''
    try:
        t = _read('README.rst')
    except UnicodeEncodeError as x:
        if t:
            i = int(str(x).split('position ')[1].split(':')[0])
            t = t[i:i+64]
        print('%s: %r' % (x, t))
        raise
    return t


def _pdf1page(pdf):
    b = _read(pdf, utf='')  # raw
    i =  b.find(b'/Type /Pages')
    if i > 0:
        i = b.find(b'/Count', i, i + 64)
        if i > 0:
            b = b[i + 6:i + 9].strip()
            return b == b'1'
    return False


def _printf(fmt, *args):
    if args:
        fmt = fmt % args
    print(fmt)


@contextmanager
def _pushd(todir):
    # <https://StackOverflow.com/questions/6194499/pushd-through-os-system>
    cwd = os.getcwd()
    os.chdir(todir)
    yield
    os.chdir(cwd)


def _read(name, utf='utf-8'):
    # read and decode all of file C{name}
    with open(name, 'rb') as f:
        t = f.read()  # all 1<<20?
    if utf and isinstance(t, bytes):
        t = t.decode(utf)
    return t


def _readlines3(name, **utf):
    t = _read(name, **utf)
    n = 0
    for s in t.split('\n'):  # os.linesep
        n += 1
        yield n, s, t


def _secs2str(secs):  # see .test/base.secs2str
    # convert secs to string
    if secs < 100:
        unit = len(_SIsecs) - 1
        while 0 < secs < 1 and unit > 0:
            secs *= 1000.0
            unit -= 1
        t = '%.3f %s' % (secs, _SIsecs[unit])
    else:
        m, s = divmod(secs, 60)
        if m < 60:
            t = '%d:%06.3f' % (int(m), s)
        else:
            h, m = divmod(int(m), 60)
            t = '%d:%02d:%06.3f' % (h, m, s)
    return t


def _selfs(py, d):
    t = _read(py)
    for c in iter('()[]{}+-=,:;'):
        t = t.replace(c, ' ')
    for w in t.split():
        if w.startswith('self.'):
            d[w] = 1 + d.get(w, 0)


def _urlower(item):
    return item[0].lower()


# <https://StackOverflow.com/questions/6194499/pushd-through-os-system>
# _os_pushd_popd_stack = list()
#
# def _os_pushdir(dirname):
#     _os_pushd_popd_stack.append(os.getcwd())
#     os.chdir(dirname)
#
# def _os_popdir():
#     os.chdir(_os_pushd_popd_stack.pop())


def _urls2(py):
    # find URLs <http...> or {http...}
    t = _read(py)
    j = 0
    while True:
        i = t.find('http', j + 1)
        if i < 0:
            break
        u, b = None, {'<': '>', '{': '}', ' ': '\n', "'": "'"}.get(t[i-1], '')
        if b:
            j = t.find(b, i + 4)
            if j > i:
                u = ''.join(t[i:j].replace('# ', '').split())
                if u.startswith('https://') or \
                   u.startswith('http://'):
                    if b == '>' and t[j+1] == '}':
                        k = t.rfind('U{', 0, i)
                        if 0 < k < i:
                            k = ' '.join(t[k+2:i-1].strip().split())
                        else:
                            k = ''
                    else:
                        k = ''
                    yield k, u
        if not u:
            _printf('bad %r url %r...', py, t[i-1:i+32])
            j = i


def _version():
    v = _version_str('pygeodesy/__init__.py')
    return '.'.join(map(str, map(int, v.split('.'))))


def _version_str(py):
    for _, s, _ in _readlines3(py):
        if s.startswith('__version__'):
            return s.split('#')[0].split('=')[-1].strip().strip("'").strip('"')
    raise ValueError('no %r in %r' % ('__version__', py))


def _versions(python):  # PYCHOK expected
    '''Invoke a different python version.
    '''
    try:
        c = python, '-m', 'test.bases'
        p = Popen(c, creationflags=0,
                   # executable   =sys.executable,
                   # shell        =True,
                     stdin        =None,
                     stdout       =PIPE,  # XXX
                     stderr       =STDOUT)  # XXX

        r = p.communicate()[0]
        if isinstance(r, bytes):  # Python 3+
            r = r.decode('utf-8')
        e = p.returncode
        if not e:
            return r.strip()
    except OSError as x:
        e = str(x)
    raise RuntimeError('%r error: %s' % (' '.join(c), e))


def _walker(start, *files):
    '''Yield all "*.py" files, deep.
    '''
    ps = ['pygeodesy', 'test'][:start]
    _join = os.path.join
    while ps:
        for p, ds, fs in os.walk(ps.pop(0)):
            # for n in sorted(ds):
            #     ps.append(_join(p, n))
            for n in sorted(fs):
                if n.endswith('.py'):
                    yield _join(p, n)
    for f in files:
        yield f


def _write(name, text, orig=''):
    if text != orig or not orig:
        with open(name, 'w') as f:
            f.write(text)


_PyGeodesy_version = 'PyGeodesy-' + _version()

if not os.getcwd().endswith('/' + _PyGeodesy_):
    sys.exit('wrong working dir')

args = sys.argv[1:] or ['check', 'docs', 'test', 'unit', 'dist', 'locs', 'tails', 'readme', 'pdf1page']
while args:  # twine or pipy # MCCABE 85
    _cmd('rm -rf .coverage')

    t = args.pop(0)
    if len(t) < 2:
        _printf('unknown target: %s', t)
        break

    if 'check'.startswith(t):
        for v in ('', '3'):
            if _cmd('pychok%s -no pygeodesy', v):
                sys.exit('pychok%s' % (v,))
            if _cmd('pychok%s -no test', v):
                sys.exit('pychok%s' % (v,))

    elif 'dist'.startswith(t) and len(t) > 1:
        _dist('pygeodesy')
        v = '3.10'  # sys.version.split('.', 1)[0]  # just '2' or '3'
        _cmd('rm -rf dist')
        # <https://docs.Python.org/2/distutils/sourcedist.html>
        # <https://docs.Python.org/3.6/distutils/sourcedist.html>
        # --formats=gztar,bztar,ztar,tar,zip --dist-dir=dist --keep-temp
        _cmd('python%s%s setup.py sdist --formats=zip,gztar,bztar', v, _OO)
        if True:  # requires the wheel package: pip install wheel
            _cmd('python%s%s -m pip install --upgrade wheel', v, _OO)
            # <https://packaging.Python.org/tutorials/distributing-packages/#pure-python-wheels>
            _cmd('python%s%s setup.py bdist_wheel --universal', v, _OO)
            _cmd('rm -rf build')
            # python3.10 -m pip wheel . -w dist  # <https://stackoverflow.com/questions/70459113/unable-to-make-a-python3-wheel-because-bdist-wheel-is-an-invalid-command>
            # python3.10 -m build --wheel  # <https://blog.ganssle.io/articles/2021/10/setup-py-deprecated.html>, <https://pypa-build.readthedocs.io/>
        _cmd('rm -rf PyGeodesy.egg-info ' + _PyGeodesy_version)

        _coverage_percentage('testcoverage/index.html', 'README.rst')

    elif 'coverage'.startswith(t) and len(t) > 1:
        # _coverage_pygeodesy_cleaner()
        _coverage_percentage('testcoverage/index.html', 'README.rst')

    elif 'docs'.startswith(t) and len(t) > 2:
        # os.environ['PYGEODESY_MAKE_DOCS'] = 'PYGEODESY_MAKE_DOCS'
        _cmd('env PYGEODESY_FOR_DOCS=1 epydoc --debug --html --no-private --no-source --name=PyGeodesy --url=https://GitHub.com/mrJean1/PyGeodesy -v pygeodesy')
        for h in iglob('html/pygeodesy*.html'):
            _html(h, *_html_from_tos)
        _html('html/epydoc.css',            ('font-size: +110%; font-style: italic;', 'font-size: +110%;'))  # non-italic signatures
        _html('html/pygeodesy-module.html', ('<b><code>Cartesian</code></b>', '<a href=pygeodesy-Cartesian-attributes-table.html><b><code>Cartesian</code></b></a>'),
                                            ('<b><code>LatLon</code></b>',    '<a href=pygeodesy-LatLon-attributes-table.html><b><code>LatLon</code></b></a>'),
                                            ('<code>height</code>', '<a href=pygeodesy.heights-module.html><code>height</code></a>'),
                                            ('<code>Geoid</code>',  '<a href=pygeodesy.geoids-module.html><code>Geoid</code></a>'))

        _cmd('rm -f html/pygeodesy-{Cartesian,LatLon}-attributes-table.html')
        _cmd('python3 make_attributes_table_html.py Cartesian >html/pygeodesy-Cartesian-attributes-table.html')
        _cmd('python3 make_attributes_table_html.py LatLon    >html/pygeodesy-LatLon-attributes-table.html')

        _cmd('rm -rf docs')  # shutil.rmtree('docs'), not os.rmdir()
        _cmd('mv html docs')

    elif 'pydocs'.startswith(t) and len(t) > 4:
        _cmd(' '.join('''env PYGEODESY_FOR_DOCS=1 pydoctor3 --project-name=PyGeodesy
                                                            --project-version=%s
                                                            --project-base-dir="."
                                                            --project-url=https://GitHub.com/mrJean1/PyGeodesy
                                                            --make-html
                                                            --html-output=pydocs
                                                            --html-viewsource-base=https://GitHub.com/mrJean1/PyGeodesy/tree/master
                                                            --intersphinx=https://docs.Python.org/3/objects.inv
                                                            --docformat=epytext ./pygeodesy'''.split()), _version())  # | grep -v Unknown | grep -v __all__
        # --intersphinx=https://docs.Python.org/3/objects.inv for external documentation links, like L{Random} in L{hausdorff}

    elif 'tables'.startswith(t) and len(t) > 4:
        _cmd('rm -f pygeodesy-{Cartesian,LatLon}-attributes.{html,pdf}')

        _cmd('python3 make_attributes_table_html.py LatLon    X >pygeodesy-LatLon-attributes.html')
        _cmd('python3 make_attributes_table_html.py Cartesian X >pygeodesy-Cartesian-attributes.html')

        _cmd('wkhtmltopdf --page-width 180.40 --page-height 530 pygeodesy-LatLon-attributes.html    pygeodesy-LatLon-attributes.pdf')  # -s A4, Letter, Legal, Tabloid, etc.
        _cmd('wkhtmltopdf --page-width 169.87 --page-height 290 pygeodesy-Cartesian-attributes.html pygeodesy-Cartesian-attributes.pdf')

        _cmd('rm -f pygeodesy-*-attributes.html')

    elif 'epytext'.startswith(t) and t > 2 and args:
        d = {}
        for f in _walker(1):  # _walker(1, 'setup.py'):
            _epytext(f, d, *args)
        print('\n'.join('%s %s' % (c, _files(fs)) for c, fs in sorted(d.items())).encode('utf-8'))
        # print('\n'.join(sorted(d.keys())).encode('utf-8'))
        break

    elif 'froms'.startswith(t) and len(t) > 1:
        _froms()

    elif 'test'.startswith(t):
        try:
            for v in ('env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=3 python3.13',
                      'env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=3 PYGEODESY_FSUM_F2PRODUCT=std PYGEODESY_GEODSOLVE=/opt/local/bin/GeodSolve PYGEODESY_INTERSECTTOOL=/opt/local/bin/IntersectTool PYGEODESY_RHUMBSOLVE=/opt/local/bin/RhumbSolve python3.12',  # geographiclib 2.0
                      'env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=0 PYGEODESY_INIT__ALL__=__all__ python3.11',  # geographiclib 2.0, numpy 1.24.2, scipy 1.10.1
                      'env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_GEOCONVERT=/opt/local/bin/GeoConvert PYGEODESY_GEODSOLVE=/opt/local/bin/GeodSolve python3.10',  # geographiclib 2.0, numpy 1.22.4, scipy 1.8.1  PYTHONDONTWRITEBYTECODE=1
#                     'env PYGEODESY_COVERAGE=0 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=0 PYTHONDONTWRITEBYTECODE=1 python3.9',  # no geographiclib, no numpy, no scipy
#                     'env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=0 PYGEODESY_FSUM_PARTIALS=0 PYGEODESY_NOTIMPLEMENTED=std python3.9',  # no geographiclib, no numpy, no scipy
#                     'env PYGEODESY_COVERAGE=0 PYGEODESY_WARNINGS=on PYGEODESY_GEOGRAPHICLIB=1 PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3.8',  # geographiclib 1.52, numpy 1.19.2, scipy 1.5.2
                      'env PYGEODESY_COVERAGE=0 PYGEODESY_WARNINGS=on PYGEODESY_GEOCONVERT=/opt/local/bin/GeoConvert PYGEODESY_GEODSOLVE=/opt/local/bin/GeodSolve PYGEODESY_INTERSECTTOOL=/opt/local/bin/IntersectTool PYGEODESY_RHUMBSOLVE=/opt/local/bin/RhumbSolve python2.7',  # geographiclib 1.50, numpy 1.16.6, scipy 1.2.2
#                     'env PYGEODESY_COVERAGE=0 PYGEODESY_WARNINGS=on PYGEODESY_LAZY_IMPORT=3 PYGEODESY_XPACKAGES=numpy,scipy,geographiclib pypy3.10',  # no geographiclib, no numpy, no scipy homebrew /opt/local/bin/pypy3.10
#                    'env PYGEODESY_COVERAGE=5 PYGEODESY_WARNINGS=on PYGEODESY_GEOGRAPHICLIB=1 PYGEODESY_LAZY_IMPORT=0 PYGEODESY_GEODSOLVE=/opt/local/bin/GeodSolve PYGEODESY_RHUMBSOLVE=/opt/local/bin/RhumbSolve PYTHONPATH=./test python3.8'
                    ):  # PYCHOK indent
                if _cmd('%s%s -W default test/run.py -results', v, _OO):
                    _printf('%s%s tests FAILED', v, _OO)
                    break
        except KeyboardInterrupt:
            pass

        v = 'python3.13'
        t = 'testcoverage (%s)' % (_PyGeodesy_version.replace('-', ' '),)
        _cmd('rm -rf testcoverage/*.*')  # *.*{css,html,js,json,png}
        _cmd('%s -m coverage html --rcfile=testcoverage.rc --title="%s"' % (v, t))
        _cmd('rm -rf .coverage')

        x = _coverage_pygeodesy_cleaner()
        if x:
            h = 'testcoverage/index__.html'
            _write(h, x)
#               assert f.name == h
            # <https://wkHTMLtoPDF.org/usage/wkhtmltopdf.txt>
            # <https://Doc.Qt.io/archives/qt-4.8/qprinter.html#PaperSize-enum>
            _cmd('wkhtmltopdf --page-width 180 --page-height 750 --disable-external-links --disable-internal-links %s testcoverage.pdf' % (h,))  # -s A4, Letter, Legal, Tabloid, etc.
            _cmd('rm -f %s' % (h,))

        _coverage_percentage('testcoverage/index.html', 'README.rst')

#       no PyPy on macOS 10.15.5 Catalina
#       for v in ('python2.7',  # 'pypy3', 'pypy2' 'intelpython3'
#                 'env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 PYGEODESY_WARNINGS=on python3.9',
#                 'env PYGEODESY_LAZY_IMPORT=0 PYGEODESY_WARNINGS=on python3.9'):
#           if _cmd('%s%s test/run.py -results', v, _OO):
#               sys.exit('%s%s tests FAILED' % (v, _OO))

    elif 'pdf1page'.startswith(t):
        if not _pdf1page('testcoverage.pdf'):
            sys.exit('%s not single page' % ('testcoverage.pdf',))

    elif 'testElevations'.startswith(t):
        for v in ('env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3',
                  'env PYGEODESY_LAZY_IMPORT=0 python3',
                # 'env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3.8',
                  'env PYGEODESY_LAZY_IMPORT=0 python3.8',
                # 'pypy3', 'pypy2',  # no PyPy on macOS 10.15.5 Catalina
                  'python2.7'):  # 'intelpython3'
            if _cmd('%s%s test/testElevations.py 5', v, _OO):
                sys.exit('%s%s tests FAILED' % (v, _OO))
            sleep(5)

    elif 'testGeodTest'.startswith(t):
        for v in ('env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3',
                  'env PYGEODESY_LAZY_IMPORT=0 python3',
                # 'env PYGEODESY_LAZY_IMPORT=0 PYTHONDONTWRITEBYTECODE=1 python3.8',
                  'env PYGEODESY_LAZY_IMPORT=3 python3.8',
                # 'pypy3', 'pypy2',  # no PyPy on macOS 10.15.5 Catalina
                  'python2.7'):  # 'intelpython3'
            if _cmd('%s%s test/testEllipsoidalGeodTest.py ../testGeodTest.dat.txt', v, _OO):
                sys.exit('%s%s tests FAILED' % (v, _OO))
            sleep(5)

    elif 'testGeoidHeights'.startswith(t):
        for v in ('env PYGEODESY_LAZY_IMPORT=3 python3',  # arm64
                # 'env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3',
                  'env PYGEODESY_LAZY_IMPORT=3 python3.8',  # arm64_x86_64
                # 'pypy3', 'pypy2',  # no PyPy on macOS 10.15.5 Catalina
                # 'python2.7',  # 'intelpython3'
                 ):  # PYCHOK indent
            for g in ('Karney', 'Pgm'):
                for e in ('egm84-15', 'egm96-5', 'egm2008-1'):
                    p = v.split()[-1].replace('p', 'P')
                    x = './testresults/testGeoid%s-%s-%s.txt' % (g, e, p)
                    # _cmd('rm -f %s', x)
                    t = time()
                    if _cmd('%s%s test/testGeoids.py -%s -dat ../testGeoids/GeoidHeights.dat -eps 1.00 ../testGeoids/%s.pgm >%s', v, _OO, g, e, x):  # '/dev/null'):
                        # sys.exit('%s%s tests FAILED' % (v, _OO))
                        _printf('%s%s tests FAILED', v, _OO)
                    t = time() - t
                    _printf('time %.3f secs for %s', t, x)

    elif 'testTMcoords'.startswith(t):
        for m in ('Exact', 'Etm', 'Utm', 'UtmUps'):
            for v in ('env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3',
                      'env PYGEODESY_LAZY_IMPORT=0 python3',
                    # 'env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3.8',
                      'env PYGEODESY_LAZY_IMPORT=0 python3.8',
                    # 'pypy3', 'pypy2',  # no PyPy on macOS 10.15.5 Catalina
                      'python2.7'):  # 'intelpython3'
                if _cmd('%s%s test/test%sTMcoords.py ../testTMcoords.dat.258K.txt', v, _OO, m):
                    sys.exit('%s%s tests FAILED' % (v, _OO))
                sleep(5)

    elif 'unit'.startswith(t):
        _gc('unit test')
        from test import TestSuite
        # <https://docs.Python.org/2/library/unittest.html>
        # <https://docs.Python.org/3.6/library/unittest.html>
        import unittest
        t = unittest.TestLoader().loadTestsFromTestCase(TestSuite)
        unittest.TextTestRunner(verbosity=2).run(t)

        for v in ('env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3',
                  'env PYGEODESY_LAZY_IMPORT=0 python3',
                  'env PYGEODESY_LAZY_IMPORT=3 PYTHONDONTWRITEBYTECODE=1 python3.8',
                  'env PYGEODESY_LAZY_IMPORT=0 python3.8',
                # 'pypy2', 'pypy3',  # no setuptools in PyPy 7.0.0, no PyPy on macOS 10.15.5 Catalina
                  'python2.7'):  # 'intelpython3'
            _cmd('%s setup.py test', v)
            _cmd('rm -rf PyGeodesy.egg-info')

    elif 'zip'.startswith(t):
        z = 'docs-%s.zip' % (_PyGeodesy_version,)
        _cmd('rm -f %s', z)
        _cmd('zip -r -9 %s index.html docs', z)

    elif t.lower() in ('twine', 'pypi'):
        if args:
            # see <https://wiki.Python.org/moin/TestPyPI>
            # and index_servers settings in file ~/.pypirc
            _PyPI, _X = 'pypitest', args.pop(0)
            if not _X.isalpha():
                sys.exit('invalid test: %s %s' % (t, _X))
        else:
            _PyPI, _X = 'pypi', ''

        _long_description()

        x = None  # only one sdist may be uploaded to PyPI
        for d in ('-py2.py3-none-any.whl', '.zip'):  # '.tar.gz',
            s = 'dist/%s%s' % (_PyGeodesy_version, d)
            d = 'dist/%s%s%s' % (_PyGeodesy_version, _X, d)
            if d != s:
                _cmd('mv %s %s', s, d)
            if os.path.exists(d):
                # _cmd('twine3 register %s -r %s' % (d, _PyPI))
                if _cmd('twine3 upload --verbose %s -r %s', d, _PyPI):
                    _printf('upload FAILED: %r', d)
            else:
                x = 'no such file: %s' % (d,)
            if s != d:
                _cmd('mv %s %s', d, s)
        if x:
            _printf(x + os.linesep)
            break

    elif t.startswith('L'):
        from pygeodesy.lazily import _all_imports
        e, p = 0, _all_imports()
#       for i, t in enumerate(sorted(p.items())):
#           print(i, t)
        for f in _walker(2, 'setup.py'):
            for n, L, t in _eLs(f, p):
                e += 1
                print('%5s %s:%s: %s ... %s' % (e, f, n, t.strip(), L))

    elif 'locs'.startswith(t):
        locs = {}
        for f in _walker(2, 'setup.py'):
            n = len(open(f).readlines())
            f = f.split(os.sep)[0]
            locs[f] = locs.get(f, 0) + n
        t = 0
        for f, n in sorted(locs.items()):
            t += n
            print('%7d %7d  %s' % (n, t, f))

    elif 'tails'.startswith(t):
        _cmd('tail -1 testresults/*.txt')

    elif 'readme'.startswith(t):
        _cmd('python3.9 -m readme_renderer README.rst -o README.html')

    elif 'urlist'.startswith(t):
        urls = {}
        n = m = 0
        for f in _walker(2, 'README.rst', 'setup.py'):
            for k, u in _urls2(f):
                if not k:
                    continue
                n += len(u)
                if k not in urls:
                    urls[k] = u
                elif urls[k] == u:
                    m += len(u)
                    continue
                t = '*' if urls[k] != u else ':'
                print('%r%s %r %s' % (k, t, u, f))
        print('duplicate %s / %s bytes or %.1f%%' % (m, n, (m * 100.0 / n)))

    elif 'urls'.startswith(t):
        try:
            from httplib import HTTPConnection as Http
            from urlparse import urlparse
        except ImportError:  # Python 3+
            from http.client import HTTPConnection as Http
            from urllib.parse import urlparse
        urls = {}
        for f in _walker(2, 'README.rst', 'setup.py'):
            for k, u in _urls2(f):
                if u in urls:
                    urls[u].append(f)
                else:
                    urls[u] = [f]
                # if k:
                #     print('%r: %r,' % (k, u))
        i = d = e = t = 0
        fs = []
        for u, f in sorted(urls.items(), key=_urlower):
            # <https://stackoverflow.com/questions/107405/
            #  how-do-you-send-a-head-http-request-in-python-2/4421712#4421712>
            p = urlparse(u.rstrip('.'))  # remove trailing '...'
            try:
                h = Http(p.netloc, timeout=4)
                h.request('HEAD', p.path)
                r = h.getresponse()
                s = r.status
                if int(s) < 400:  # 301 Moved,
                    s = 'OK'
                else:
                    s = '%s *****' % (s,)
                    e += 1
            except Exception as x:
                s = '%r *****' % (x,)
                e += 1
                if s.startswith('timeout'):
                    t += 1
            i += 1
            d += len(f)
            _printf('%5d: <%s> %s %s', i, u, s, _files(f))
            if s != 'OK':
                fs.append(i)
        d -= i  # duplicates
        if e:
            f = 'FAILED'
            if t:
                f += ', incl. %d timed out' % (t,)
            if fs:
                f = '%s (%s)' % (f, ', '.join(map(str, fs)))
            s = 's' if e > 1 else ''
        else:
            f = 'passed'
            e = 'all'
            s = 's'
        _printf('\n%s urls check%s %s, %d duplicates', e, s, f, d)

    elif t in ('xtgz',) and args:
        p = _PyGeodesy_
        v = p + args.pop(0)
        with _pushd('./..'):
            _cmd('env COPYFILE_DISABLE=true COPY_EXTENDED_ATTRIBUTES_DISABLE=true tar -cvf %s %s', v, p)
            _cmd('gzip -9 -S .tgz %s', v)
            _cmd('ls -l *.tgz')
        break

    elif 'exported'.startswith(t):
        # print all exported, public classes, constants, functions and modules
        import pygeodesy as P
        from inspect import isfunction, ismodule

        _Enum = P.named._NamedEnum

        base = P.__name__
        dups = 0
        line = 0
        prev = ''

        def _export(mod, nam, val=''):
            global dups, line, prev
            nam = mod + '.' + nam
            if nam == prev:
                dups += 1
            else:
                prev = nam
            nam = nam + ' ' + val
            if len(nam) > 160:
                nam = nam[:80] + '....' + nam[-80:]
            line += 1
            _printf('%3d %s', line, nam)

        def _repr(o, n):
            if isfunction(o):
                return "<function '%s.%s'>" % (o.__module__, n)
            elif ismodule(o):
                return '<module %r %s>' % (o.__name__, o.__version__)
            elif isinstance(o, str):
                return "%s('%s')" % (type(o).__name__, o)
            else:
                r = repr(o)
                if not (r.startswith('<') and r.endswith('>')):
                    r = '%s(%s)' % (type(o).__name__, r)
                return r

        for n in sorted(P.__all__, key=str.lower):
            o = getattr(P, n)
            if isinstance(o, _Enum):
                for n in repr(o).split(',\n'):
                    _export(base, n)
            elif ismodule(o):
                _export(base, n, '<module %r>' % (o.__name__,))
                if o in (P.ellipsoidalKarney, P.ellipsoidalNvector, P.ellipsoidalVincenty,
                         P.sphericalNvector, P.sphericalTrigonometry):
                    for n in sorted(o.__all__, key=str.lower):
                        r = _repr(getattr(o, n), n)
                        _export(o.__name__, n, r)
            else:
                r = _repr(o, n).replace(' ' + n + ' ', ' ') \
                               .replace(" '" + n + "' ", ' ')
                _export(base, n, r)
        _printf('--- %s %s (%s duplicates)', P.__name__, P.version, dups or 'no')

    elif 'self'.startswith(t):
        n, d = 0, {}
        while args:
            _selfs(args.pop(0), d)
            n += 1
        for t in sorted(d.items()):
            print(' %s: %sX' % t)
        print('\n %d selfs, %dX, %d files' % (len(d), sum(d.values()), n))

    elif 'versions'.startswith(t) and args:
        _printf(os.linesep + _versions(args.pop(0)))

    elif '__version__'.startswith(t) and len(t) > 2:
        import re
        _valid = '[0-9][0-9][.][0-9][0-9][.][0-9][0-9]'
        _valid =  re.compile(_valid).match

        for f in _walker(2, 'setup.py', 'make'):
            v = _version_str(f)
            b = '' if len(v) == 8 and _valid(v) else 'INVALID!'
            print('%s = %r in %s %s' % ('__version__', v, f, b))
            if b:
                break

    elif '__all__'.startswith(t) and len(t) > 4:
        _all_locals(*args)

    else:
        _printf('unknown target: %s', t)
        break

_gc('make, Python ' + sys.version.split()[0])
