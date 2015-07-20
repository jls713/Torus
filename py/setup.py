#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

# os.environ["CC"] = "g++-4.9"
# os.environ["CXX"] = "g++-4.9"

TorusPath = '../'

setup(name="Torus_py",
      ext_modules=[
          Extension("Torus_py", ["src/Torus_py.cpp"],
                    include_dirs=[TorusPath + 'src',
                                  TorusPath + 'src/utils',
                                  np.get_include()],
                    library_dirs=[TorusPath+'obj/'],
                    libraries=[
                        'Torus', 'Other', 'Pot', 'WD', 'ebf_cpp',
                        'boost_python', 'm', 'c++'],
                    extra_compile_args=['-fPIC', '-shared', '-std=c++0x',
                                        '-Wall', '-O3', '-Ofast',
                                        '-ffast-math', '-fPIC', '-fopenmp'])
      ])
