from distutils.core import setup, Extension
import numpy as np
setup(name='summing_numpy', version='1.0', include_dirs = [np.get_include()], ext_modules=[Extension('summing_numpy', sources = ['/home/evgenii/crossing/crossing.c'])])
