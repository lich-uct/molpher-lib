
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
import shutil

from jinja2 import Template

# The short X.Y version.
current_version = 'v' + os.environ['VERSION']
versions = [ x for x in os.listdir(".") if x.startswith("v")]

if current_version in versions:
    print("WARNING: Docs for this version already built. Rewriting...")
    shutil.rmtree(current_version)
else:
    versions.append(current_version)

INDEX_TEMPLATE_PATH = "index.html.template"
INDEX_PATH = "index.html"

versions.sort()

template = None
with open(INDEX_TEMPLATE_PATH, "r", encoding="utf-8") as template_file:
    template = Template(template_file.read())

with open(INDEX_PATH, "w", encoding='utf-8') as index_file:
    index_file.write(template.render(
        versions=versions
    ))