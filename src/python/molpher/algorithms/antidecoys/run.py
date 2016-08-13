import os
import pickle

from molpher.algorithms.commons import timeit
from .utils import compute_anti_fp
from .pathfinder import PathFinder

def run(settings, paths_count):
    storage_dir = settings.storage_dir
    antifp_path = os.path.join(storage_dir, 'antifingerprint.pickle')
    paths_path = os.path.join(storage_dir, 'paths.pickle')

    paths = []
    antifp = None
    if os.path.exists(paths_path):
        pickled_file = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_file))
        pickled_file.close()
    if os.path.exists(antifp_path):
        pickled_file = open(antifp_path, mode='rb')
        antifp = pickle.load(pickled_file)
        pickled_file.close()

    pathfinder = PathFinder(
        settings
        , antifingerprint=antifp
    )
    for i in range(paths_count):
        # find a path
        exec_time = timeit(pathfinder)
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(exec_time, i + 1))
        print('----------------------------------------------------------------')

        if pathfinder.path:
            antifp = compute_anti_fp(pathfinder.path, settings.signature_factory, antifp)
        pathfinder = PathFinder(
            settings
            , antifingerprint=antifp
        )

        # pickle the results for future use
        if antifp:
            pickled_file = open(antifp_path, mode='wb')
            pickle.dump(antifp, pickled_file)
            pickled_file.close()
        pickled_file = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_file)
        pickled_file.close()