from setuptools import setup, Extension, find_packages
import os
from version import VERSION

my_cpp_module = Extension('molpher.swig_wrappers._core',
                           sources=['src/python/molpher_wrap.cpp'],
                           include_dirs = [
                                os.path.abspath('include/')
                                , os.path.abspath('../backend/')
                                , os.path.abspath('../common/')
                                , os.path.abspath('../dependencies/rdkit/Code')
                                , os.path.abspath('../dependencies/boost')
                                , os.path.abspath('../dependencies/tbb/include')
                                ],
                           library_dirs=[os.path.join(os.getcwd(), 'lib/')],
                           libraries=['molpher'],
                           runtime_library_dirs=[os.path.join(os.getcwd(), 'lib/')],
                           extra_compile_args=['-std=c++11'],
                           language='c++'
                           )

setup (name = 'molpher',
       version = VERSION,
       author      = "Martin Sicho",
       description = "Molpher-lib Python wrappers and library.",
       package_dir={'': 'python/'},
       packages=find_packages('python/'),
       ext_modules = [my_cpp_module,],
       package_data = {
            #'molpher.swig_wrappers': ['lib/*.so', 'lib/*.so.2', '*.dat'],
            'molpher.swig_wrappers': ['*.dat'],
            'molpher.swig_wrappers.tests' : ['test_files/*']
        },
       test_suite="tests"
       )
