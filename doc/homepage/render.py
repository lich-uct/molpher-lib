
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

from jinja2 import Template, Environment, FileSystemLoader

from pkg_resources import parse_version

def render_template(temp, out_path, ctx):
    with open(out_path, "w", encoding='utf-8') as index_file:
        index_file.write(temp.render(
            **ctx
        ))

def prepare_versions():
    current_version = 'v' + os.environ['VERSION']
    versions = [ x for x in os.listdir(".") if x.startswith("v")]
    if current_version in versions:
        print("WARNING: Docs for this version already built. Rewriting...")
        shutil.rmtree(current_version)
    else:
        versions.append(current_version)
    versions.sort(key=parse_version) # sort from newest to oldest
    versions.reverse()
    newest = None
    for item in versions:
        if not 'dev' in item:
            newest = item
            break
    return {
        'versions' : versions
        , 'newest' : newest
    }

# initialize context with version info
context = prepare_versions()

# prepare environment
env = Environment()
env.loader = FileSystemLoader('.')

# render the homepage
context['active'] = 'index'
render_template(env.get_template('index.html.template'), "index.html", context)

# render the examples page
context['active'] = 'examples'
render_template(env.get_template('examples.html.template'), "examples.html", context)