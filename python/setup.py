from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ( "pycube",
              ['pycube.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args=["-std=c++17"],
              libraries = ['cube'],
              include_dirs = ['../cube', '../lib'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pycube',
    ext_modules = cythonize(extensions, language_level = "3"),
    )
