
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
import sys
from pkg_resources import resource_filename

PYTHON_VERSION_MAJOR = sys.version_info.major
PYTHON_VERSION_MINOR = sys.version_info.minor
PYTHON_VERSION = '{0}.{1}'.format(PYTHON_VERSION_MAJOR, PYTHON_VERSION_MINOR)

def get_info():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../'))
    version_path = os.path.join(base_dir, "VERSION.TXT")
    build_path = os.path.join(base_dir, "BUILD.TXT")
    if not os.path.exists(version_path):
        version_path = resource_filename('molpher', 'VERSION.TXT')
        build_path = resource_filename('molpher', 'BUILD.TXT')

    version = None
    with open(version_path, 'r') as versionfile:
        version = versionfile.read().strip()

    build = None
    with open(build_path, 'r') as buildfile:
        build = buildfile.read().strip()

    return version, build

VERSION, BUILD_NUMBER = get_info()