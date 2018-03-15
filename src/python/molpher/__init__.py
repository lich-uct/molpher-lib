"""
When this module is imported,
some actions are automatically
performed (see below).

"""

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

from pkg_resources import resource_filename

import molpher.swig_wrappers.core as wrappers

import os

_SHARED_FOLDER = os.path.join(
   os.path.abspath(os.path.dirname(__file__)) # this module directory
   , '../../../../../share/molpher-lib' # path to the share directory with resources
)

_SA_SEARCH_PATHS = [
   os.path.join(_SHARED_FOLDER, 'SAScore.dat')
   , resource_filename('molpher.swig_wrappers', 'SAScore.dat')
   , os.path.join(os.path.abspath(os.path.dirname(__file__)), '../../../res/SAScore.dat')
]

def load_SAScore(path):
   """
   Loads the data file used for the computation of the syntetic feasibility scores.
   This is performed automatically when the :py:mod:`molpher` package is imported.
   Use it only if you want to use a data file different from the dafault one.

   :param path: path to the `SAScore.dat` file
   :type path: :py:class:`str`
   :return: :py:obj:`None`
   """

   # print("Loading data from:", path)
   if not os.path.exists(path):
      raise Exception("No SAScore data file not found at location:", path)
   wrappers.load_data_from(path)
   return True

# print("Initializing Molpher-lib...")
file_found = False
for path in _SA_SEARCH_PATHS:
   try:
      file_found = load_SAScore(path)
      break
   except Exception as exp:
      continue

if not file_found: raise Exception("No SAScore data file found. Locations tried:", _SA_SEARCH_PATHS)
# print("Molpher-lib initialized.")