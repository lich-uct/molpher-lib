from setuptools import setup, Extension, find_packages
import os
from version import VERSION

my_cpp_module = Extension('molpher.swig_wrappers._core',
                           sources=['src/swig/molpher_wrap.cpp'],
                           include_dirs = [os.path.abspath('include/')],
                           library_dirs=[os.path.abspath('dist/lib/')],
                           libraries=['molpher'],
                           runtime_library_dirs=[os.path.abspath('dist/lib/')],
                           extra_compile_args=['-std=c++11'],
                           language='c++'
                           )

setup (name = 'molpher',
       version = VERSION,
       author      = "Martin Sicho",
       description = "Molpher-lib Python wrappers and library.",
       package_dir={'': 'src/python/'},
       packages=find_packages('src/python/'),
       ext_modules = [my_cpp_module,],
       package_data = {
            #'molpher.swig_wrappers': ['lib/*.so', 'lib/*.so.2', '*.dat'],
            'molpher.swig_wrappers': ['*.dat'],
            'molpher.swig_wrappers.tests' : ['test_files/*'],
            'molpher.examples' : ['*.xml']
        }
       #test_suite="tests"
       )
