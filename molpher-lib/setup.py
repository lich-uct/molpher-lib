from setuptools import setup, Extension


my_cpp_module = Extension('molpher._core',
                           sources=['src/python/molpher_wrap.cpp'],
                           include_dirs = [
                                'include/'
                                , '../backend/'
                                , '../common/'
                                , '../dependencies/rdkit/Code'
                                , '../dependencies/boost'
                                , '../dependencies/tbb/include'
                                ],
                           library_dirs=['lib/'],
                           libraries=['molpher'],
                           runtime_library_dirs=['lib/'],
                           extra_compile_args=['-std=c++11'],
                           language='c++'
                           )

setup (name = 'molpher',
       version = '0.1',
       author      = "Martin Sicho",
       description = "Molpher Python wrappers.",
       package_dir={'molpher': 'python/molpher'},
       packages=['molpher'],
       ext_modules = [my_cpp_module,],
       )
