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
TARGETS = ["molpher-lib"]
CONDA_FLAGS = {
    "molpher-lib" : "-c rdkit --croot /tmp/conda-bld/"
}
PYTHON_VERSIONS = ['2.6', '2.7',  '3.3',  '3.4', '3.5', '3.6', '3.7'] if BUILD_ALL_PYTHON else ['3']
VERSION = version.VERSION
BUILD_NUMBER = version.BUILD_NUMBER
LICENSE_FILE_NAME = "LICENSE"

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
            ))

        os.chdir(os.path.join(PACKAGE_DIR, "../"))
        os.environ['BASE_DIR'] = BASE_DIR
        print(os.environ['PATH'])
        copyfile(os.path.join(BASE_DIR, "LICENSE.md"), os.path.join(PACKAGE_DIR, LICENSE_FILE_NAME))
        ret = subprocess.call("conda build {0} --python {1} {2}".format(target, python_version, CONDA_FLAGS[target]), env=os.environ, shell=True)

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