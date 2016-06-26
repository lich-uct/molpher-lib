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
from jinja2 import Template
import imp
version = imp.load_source('module.name', '../version.py')

BUILD_ALL_PYTHON = False
TARGETS = ["tbb", "molpher-lib"]
PYTHON_VERSIONS = ['2.6', '2.7',  '3.3',  '3.4', '3.5'] if BUILD_ALL_PYTHON else ['3.5', '2.7']
VERSION = version.VERSION
BUILD_NUMBER = version.BUILD_NUMBER
LICENSE_FILE_NAME = "LICENSE"
TBB_LICENSE_FILE = "deps/tbb/COPYING"
TBB_VERSION = "4.2"

for target in TARGETS:
    BASE_DIR = os.path.dirname(os.path.abspath(os.path.join(__file__, '..')))
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
                , version=VERSION + '_r' + BUILD_NUMBER if BUILD_NUMBER != '0' else VERSION
                , build_number=BUILD_NUMBER
                , license_file=LICENSE_FILE_NAME
                , build_string="py{0}_{1}".format(python_version.replace('.', ''), BUILD_NUMBER)
                , python_spec="python {0}*".format(python_version)
                , tbb_version=TBB_VERSION
                , tbb_spec = "tbb {0}*".format(TBB_VERSION)
            ))

        os.chdir(os.path.join(PACKAGE_DIR, "../"))
        os.environ['BASE_DIR'] = BASE_DIR
        print(os.environ['PATH'])
        if target == 'tbb':
            # do not set Python version for tbb
            copyfile(os.path.join(BASE_DIR, TBB_LICENSE_FILE), os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
            ret = subprocess.call("conda build {0}".format(target), env=os.environ, shell=True)
        else:
            copyfile(os.path.join(BASE_DIR, "LICENSE.md"), os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
            ret = subprocess.call("conda build {0} --python {1}".format(target, python_version), env=os.environ, shell=True)

        # return non-zero if something went wrong
        if ret != 0:
            print("An error has occured during conda build...")
            exit(1)

        # remove the generated and copied files and go back to BASE_DIR
        os.remove(METAFILE_PATH)
        os.remove(os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
        os.chdir(BASE_DIR)

        # make sure that the tbb package is only built once
        if target == 'tbb':
            break
