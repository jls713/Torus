#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

os.environ["CC"] = "/opt/ioa/software/gcc/4.7.2/bin/g++"
os.environ["CXX"] = "/opt/ioa/software/gcc/4.7.2/bin/g++"
BOOSTINCPATH = "/opt/ioa/software/boost/1.55.0/include"
BOOSTLIBPATH = "/opt/ioa/software/boost/1.55.0/lib"
EBFINC    = "/data/jls/libebf_c_cpp-0.0.3/include/"
EBFLIB    = "/data/jls/libebf_c_cpp-0.0.3/lib/"
PYTHONINCPATH= "/opt/ioa/software/python/2.7.8/include/python2.7/"
PYTHONLIBPATH= "/opt/ioa/software/python/2.7.8/lib/"

TorusPath = '../'

setup(name="Torus_py",
      ext_modules=[
          Extension("Torus_py", ["src/Torus_py.cpp"],
                    include_dirs=[TorusPath + 'src',
                                  TorusPath + 'src/utils',
                                  np.get_include(),BOOSTINCPATH,EBFINC,PYTHONINCPATH],
                    library_dirs=[TorusPath+'obj/',BOOSTLIBPATH,EBFLIB,PYTHONLIBPATH],
                    libraries=[
                        'Torus', 'Other', 'Pot', 'WD', 'ebf_cpp',
                        'boost_python', 'm'],
                    extra_compile_args=['-fPIC', '-shared', '-std=c++0x',
                                        '-Wall', '-O3',
                                        '-ffast-math','-fPIC', '-fopenmp'])
      ])
