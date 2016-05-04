from pkg_resources import resource_filename

import molpher.swig_wrappers.core as wrappers

def load_SAScore(path):
   """
   Loads the data file used in computation of the syntetic feasibility scores.
   This is performed automatically when the :py:mod:`molpher` package is imported.
   Use it only if you want to use a data file different from the dafault one.

   :param path: path to the `SAScore.dat` file
   :type path: :py:class:`str`
   :return: :py:obj:`None`
   """

   print("Loading data from:", path)
   wrappers.load_data_from(path)

load_SAScore(resource_filename('molpher.swig_wrappers', 'SAScore.dat'))