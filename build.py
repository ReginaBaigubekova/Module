from distutils.core import setup, Extension

setup(name='SLEsolver', version='1.0', ext_modules=[Extension('SLEsolver', ['moduleSLE.c'])])