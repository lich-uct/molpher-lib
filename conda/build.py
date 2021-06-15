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
import subprocess
from shutil import copyfile

from jinja2 import Template
from importlib.machinery import SourceFileLoader

JOBS = 6
BASE_DIR = os.path.dirname(os.path.abspath(os.path.join(__file__, '..')))
pythons = [f'3.{x}' for x in range(6,10)]

version_module = SourceFileLoader('version', os.path.join(BASE_DIR, 'src/python/molpher/version.py')).load_module()
version = version_module.VERSION
build = version_module.BUILD_NUMBER

def generate_config_file(template, output, python_version):
    with open(template, "r", encoding="utf-8") as metafile_template:
        print('Using template:', os.path.abspath(template))
        template = Template(metafile_template.read())
    with open(output, "w", encoding='utf-8') as metafile:
        print('Generating configuration file:', os.path.abspath(output))
        metafile.write(template.render(
            version=version
            , build_number=build
            , python_version=python_version
        ))
        print('Done.')

for python_version in pythons:
    print(f'Preparing Python version: {python_version}')
    generate_config_file(
        "./molpher-lib/meta.yaml.template",
        "./molpher-lib/meta.yaml",
        python_version=python_version
    )
    os.environ['PYTHON_VERSION'] = python_version
    os.environ['JOBS'] = str(JOBS)
    os.environ['BASE_DIR'] = str(BASE_DIR)
    copyfile(os.path.join(BASE_DIR, "LICENSE.md"), os.path.join('./molpher-lib/', "LICENSE.md"))

    ret = subprocess.call("./build.sh", env=os.environ, shell=True)
    if ret != 0:
        print("An error has occured during conda build...")
        exit(1)

    print('Done.')