import os
import pickle

from .settings import AntidecoysSettings
from .utils import timeit, compute_anti_fp
from .pathfinders import PathFinder

def run(
        source
        , target
        , verbose=True
        , settings=AntidecoysSettings()
):
    storage_dir = settings.storage_dir
    if not os.path.exists(storage_dir):
        os.mkdir(storage_dir)
    antifp_path = os.path.join(storage_dir, 'antifingerprint.pickle')
    paths_path = os.path.join(storage_dir, 'paths.pickle')

    paths = []
    paths_antifp = None
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
        source
        , target
        , verbose=verbose
        , antifingerprint=antifp
        , settings=settings
    )
    for i in range(settings.max_paths):
        # find a path
        exec_time = timeit(pathfinder)
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(exec_time, i + 1))
        print('----------------------------------------------------------------')

        if pathfinder.path:
            antifp = compute_anti_fp(pathfinder.path, settings.signature_factory, antifp)
        pathfinder = PathFinder(
            source
            , target
            , verbose=verbose
            , settings=settings
            , antifingerprint=antifp)

        # pickle the results for future use
        if antifp:
            pickled_file = open(antifp_path, mode='wb')
            pickle.dump(antifp, pickled_file)
            pickled_file.close()
        pickled_file = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_file)
        pickled_file.close()