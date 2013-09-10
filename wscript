APPNAME = 'atm'
VERSION = '1.0'

top = '.'
out = 'build'

def options(opt):
    opt.load('compiler_cxx boost')

def configure(conf):
    conf.load('compiler_cxx boost')
    conf.check_boost()

    # Boost
    conf.check_boost(lib='timer', mandatory=True)

    conf.define('APP_NAME', APPNAME)
    conf.define('APP_VERSION', VERSION)
    conf.write_config_header('config.h')

    # base environment
    base_env = conf.env

    # an environment for "debug" variant
    conf.setenv('debug', base_env)
    conf.env.CFLAGS   = [              '-Wall', '-g']
    conf.env.CXXFLAGS = ['-std=c++11', '-Wall', '-g']

    # an environment for "release" variant
    conf.setenv('release', base_env)
    conf.env.CFLAGS   = [              '-Wall', '-march=native', '-O3', '-DNDEBUG']
    conf.env.CXXFLAGS = ['-std=c++11', '-Wall', '-march=native', '-O3', '-DNDEBUG']

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
