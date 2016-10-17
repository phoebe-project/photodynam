from numpy.distutils.core import setup, Extension

# Set to true if you want to link against electric fence:
CDEBUG = False

libraries = []
if CDEBUG:
    libraries += ['efence']

ext_modules = [
    Extension('photodynam',
      sources = ['./pywrapper/libphotodynam.cpp'],
      extra_compile_args = ["-std=c++11", "-I./include"]),
]

setup (name = 'photodynam',
       version = 'devel',
       description = "python wrapper of a modified version of Josh Carter's photodynam",
       packages = [],
       install_requires=['numpy'],
       ext_modules = ext_modules)