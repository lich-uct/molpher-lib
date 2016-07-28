import os
import pickle

from .settings import STORAGE_DIR, PATHS_TO_FIND
from .utils import timeit, antifingerprint_from_paths


def run(source, target, antidecoys_filter=None, use_path_antifp=True, verbose=True):
    antifp_path = os.path.join(STORAGE_DIR, 'antifingerprint.pickle')
    paths_path = os.path.join(STORAGE_DIR, 'paths.pickle')
    from .pathfinders import BidirectionalPathFinder

    paths = []
    path_antifp = None
    if os.path.exists(paths_path):
        pickled_paths = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_paths))
        pickled_paths.close()
    if os.path.exists(antifp_path): # FIXME: path_antifp and antifp should be serialized to different files
        pickled_antifp = open(antifp_path, mode='rb')
        if antidecoys_filter:
            antidecoys_filter.antifingerprint = pickle.load(pickled_antifp)
        if use_path_antifp:
            path_antifp = pickle.load(pickled_antifp)
        pickled_antifp.close()

    pathfinder = BidirectionalPathFinder(source, target, verbose=verbose, path_antifingerprint=path_antifp)
    for i in range(PATHS_TO_FIND):
        # find a path
        exec_time = timeit(pathfinder)
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(exec_time, i + 1))

        # reset the pathfinder for another search
        if pathfinder.path:
            # compute and save new antifingerprint
            if antidecoys_filter:
                antidecoys_filter.antifingerprint = antifingerprint_from_paths(pathfinder.path, antidecoys_filter.antifingerprint)
            if use_path_antifp:
                path_antifp = antifingerprint_from_paths(pathfinder.path, path_antifp)
                pathfinder.path_antifingerprint = path_antifp
            pathfinder.reset()
        else:
            pathfinder = BidirectionalPathFinder(source, target, verbose=verbose, path_antifingerprint=path_antifp)

        # pickle the results for future use
        pickled_antifp = open(antifp_path, mode='wb') # FIXME: path_antifp and antifp should be serialized to different files
        if antidecoys_filter:
            pickle.dump(antidecoys_filter.antifingerprint, pickled_antifp)
        if use_path_antifp:
            pickle.dump(path_antifp, pickled_antifp)
        pickled_antifp.close()
        pickled_paths = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_paths)
        pickled_paths.close()