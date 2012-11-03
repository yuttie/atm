APPNAME = 'p4'
VERSION = '0.6.2'

top = '.'
out = 'build'

def options(opt):
    opt.load('compiler_cxx boost')

def configure(conf):
    conf.load('compiler_cxx boost')
    conf.check_boost()

    # base environment
    base_env = conf.env

    # an environment for "debug" variant
    conf.setenv('debug', base_env)
    conf.env.CFLAGS   = ['-g', '-ggdb', '-Wall']
    conf.env.CXXFLAGS = ['-g', '-ggdb', '-Wall', '-std=c++0x']

    # an environment for "release" variant
    conf.setenv('release', base_env)
    conf.env.CFLAGS   = ['-march=native', '-O2', '-DNDEBUG', '-Wall']
    conf.env.CXXFLAGS = ['-march=native', '-O2', '-DNDEBUG', '-Wall', '-std=c++0x']

def build(bld):
    bld.recurse('src')

# define "debug" and "release" variants
from waflib.Build import BuildContext, CleanContext, \
        InstallContext, UninstallContext

for c in (BuildContext, CleanContext, InstallContext, UninstallContext):
    name = c.__name__.replace('Context','').lower()
    # declare a command for "debug" variant
    class tmp(c):
        cmd = name
        variant = 'debug'
    # declare a command for "release" variant
    class tmp(c):
        cmd = name + '_release'
        variant = 'release'
