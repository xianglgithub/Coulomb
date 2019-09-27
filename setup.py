
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys
import os
import shutil


# clean previous build
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if (name.startswith("libCoulombCH3I") and not(name.endswith(".pyx") or name.endswith(".pxd") or name.endswith(".a"))):
            os.remove(os.path.join(root, name))
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)

ext_modules = [Extension("libCoulombCH3I",
                     ["libCoulombCH3I.pyx"],
                     libraries=["gsl","gslcblas","m"],
                     language='c++',                         
                     extra_objects=["libCoulombCH3I.a"],   
                     extra_compile_args=["-I/usr/includes/"],
                     extra_link_args=["-L/usr/lib64"]
                     )]

setup(
  name = 'libCoulombCH3I',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)





