
# Copyright (c) 2016 Martin Sicho
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, Extension, find_packages
import os
from version import VERSION, BUILD_NUMBER

prefix_include_dir = os.path.join(os.environ.get('PREFIX'), 'include/molpher-lib/') if os.environ.get('PREFIX') else None
conda_prefix_include_dir = os.path.join(os.environ.get('CONDA_PREFIX'), 'include/molpher-lib/') if os.environ.get('CONDA_PREFIX') else None

include_dirs = [
    os.path.abspath('dist/include/molpher-lib/')
]
if prefix_include_dir:
    include_dirs += [prefix_include_dir]
if conda_prefix_include_dir:
    include_dirs += [conda_prefix_include_dir]

molpher_cpp_module = Extension('molpher.swig_wrappers._core',
                           sources=[os.path.abspath('src/swig/molpher_wrap.cpp')],
                           include_dirs = include_dirs,
                           library_dirs = os.environ["LD_LIBRARY_PATH"].split(":") if "LD_LIBRARY_PATH" in os.environ else None,
                           libraries=['molpher'],
                           extra_compile_args=['-std=gnu++1z'],
                           language='c++'
                           )

setup (name = 'molpher',
       version = VERSION + '-r' + BUILD_NUMBER if BUILD_NUMBER != '0' else VERSION,
       author      = "Martin Sicho",
       description = "Molpher-lib Python wrappers and library.",
       package_dir={'': 'src/python/'},
       packages=find_packages('src/python/'),
       ext_modules = [molpher_cpp_module,],
       package_data = {
            'molpher.core.tests' : ['test_files/*'],
            'molpher.examples' : ['*.xml', '*.sdf'],
            'molpher' : ['SAScore.dat']
        },
       test_suite="molpher.core.tests"
       )
