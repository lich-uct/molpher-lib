import os
import pickle

from molpher.algorithms.classic.pathfinder import ClassicPathFinder


def run(settings, iter_callback=None):
    """
    Run a search with the given settings and return the path found.

    :param iter_callback: callable that is run after each iteration, the exploration tree is supplied as its sole parameter
    :param settings: `Settings` instance
    :return: path as a `list` of `str`
    """

    pathfinder = ClassicPathFinder(settings, iter_callabck=iter_callback)
    pathfinder()

    # pickle the resulting path
    pickled_path = open(os.path.join(settings.storage_dir, 'path.pickle'), mode='wb')
    pickle.dump(pathfinder.path, pickled_path)
    pickled_path.close()

    return pathfinder.path