import os
import pickle

from molpher.algorithms.functions import timeit
from .utils import compute_anti_fp
from .pathfinder import AntidecoysPathFinder

def run(settings, paths_count):
    """

    :param settings: instance of `AntidecoysSettings` specifying search parameters
    :param paths_count: number of paths to generate
    :return: a `list` of lists of SMILES strings (each path is represented by a list of SMILES strings)
    """

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

    pathfinder = AntidecoysPathFinder(
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
        pathfinder = AntidecoysPathFinder(
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

    return paths