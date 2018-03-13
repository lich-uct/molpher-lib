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

import os
from shutil import copyfile
import subprocess

import shutil
from jinja2 import Template
import imp

BASE_DIR = os.path.dirname(os.path.abspath(os.path.join(__file__, '..')))
version = imp.load_source('module.name', os.path.join(BASE_DIR, 'version.py'))
BUILD_ALL_PYTHON = False
# TARGETS = ["molpher_tbb", "molpher_boost", "molpher_rdkit", "molpher-lib"]
TARGETS = ["molpher-lib"]
PYTHON_VERSIONS = ['2.6', '2.7',  '3.3',  '3.4', '3.5', '3.6'] if BUILD_ALL_PYTHON else ['3', '2']
VERSION = version.VERSION
BUILD_NUMBER = version.BUILD_NUMBER
LICENSE_FILE_NAME = "LICENSE"
# TBB_LICENSE_FILE = "deps/tbb/COPYING"
# BOOST_LICENSE_FILE = "deps/boost/LICENSE_1_0.txt"
# RDKIT_LICENSE_FILE = "deps/rdkit/license.txt"
# TBB_VERSION = "4.2"
# BOOST_VERSION = "1.63"
# RDKIT_VERSION = "2017.09.1"

# LICENSES_DICT = {
#     "molpher_tbb" : TBB_LICENSE_FILE
#     , "molpher_boost" : BOOST_LICENSE_FILE
#     , "molpher_rdkit" : RDKIT_LICENSE_FILE
# }

# remove previously generated files
shutil.rmtree(os.path.join(BASE_DIR, 'build'), ignore_errors=True)
shutil.rmtree(os.path.join(BASE_DIR, 'dist'), ignore_errors=True)

for target in TARGETS:
    os.chdir(BASE_DIR)
    for python_version in PYTHON_VERSIONS:
        PACKAGE_DIR = os.path.join(BASE_DIR, "conda/{0}/{0}/".format(target))
        METAFILE_TEMPLATE_PATH = os.path.join(PACKAGE_DIR, "../meta.yaml.template")
        METAFILE_PATH = os.path.join(PACKAGE_DIR, "meta.yaml")

        template = None
        with open(METAFILE_TEMPLATE_PATH, "r", encoding="utf-8") as metafile_template:
            template = Template(metafile_template.read())

        with open(METAFILE_PATH, "w", encoding='utf-8') as metafile:
            metafile.write(template.render(
                target=target
                , version=VERSION
                , build_number=BUILD_NUMBER
                , license_file=LICENSE_FILE_NAME
                , build_string="py{0}_{1}".format(python_version.replace('.', ''), BUILD_NUMBER)
                , python_spec="python {0}*".format(python_version)
                # , tbb_version=TBB_VERSION
                # , tbb_spec = "{1} {0}".format(TBB_VERSION, TARGETS[0])
                # , boost_version=BOOST_VERSION
                # , boost_spec = "{1} {0}*".format(BOOST_VERSION, TARGETS[1])
                # , rdkit_version=RDKIT_VERSION
                # , rdkit_spec = "{1} {0}".format(RDKIT_VERSION, TARGETS[2])
            ))

        os.chdir(os.path.join(PACKAGE_DIR, "../"))
        os.environ['BASE_DIR'] = BASE_DIR
        print(os.environ['PATH'])
        if target in TARGETS[:-1]:
            # do not set Python version for tbb
            # copyfile(os.path.join(BASE_DIR, LICENSES_DICT[target]), os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
            ret = subprocess.call("conda build {0} --croot /tmp/conda-bld/".format(target), env=os.environ, shell=True)
        else:
            copyfile(os.path.join(BASE_DIR, "LICENSE.md"), os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
            ret = subprocess.call("conda build -c lich -c conda-forge -c rdkit {0} --python {1} --croot /tmp/conda-bld/".format(target, python_version), env=os.environ, shell=True)

        # return non-zero if something went wrong
        if ret != 0:
            print("An error has occured during conda build...")
            exit(1)

        # remove the generated and copied files and go back to BASE_DIR
        os.remove(METAFILE_PATH)
        os.remove(os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
        os.chdir(BASE_DIR)

        # the dependencies are only built once
        if target in TARGETS[:-1]:
            break