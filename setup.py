from setuptools import Extension, setup
import numpy

# Set to true if you want to link against electric fence:
CDEBUG = False

libraries = []
if CDEBUG:
    libraries += ['efence']

ext_modules = [
    Extension(
        'photodynam',
        sources=[
            './pywrapper/libphotodynam.cpp'
        ],
        extra_compile_args=["-std=c++11", "-I./include"],
        include_dirs=[numpy.get_include()]
    ),
]

setup(ext_modules=ext_modules)
