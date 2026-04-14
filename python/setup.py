from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ( "pyatom",
              ['pyatom.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args=["-std=c++11"],
              libraries = ['atom'],
              include_dirs = ['../atmosphere', '../hydrosphere', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyatom',
    ext_modules = cythonize(extensions, language_level = "2"),
    )
