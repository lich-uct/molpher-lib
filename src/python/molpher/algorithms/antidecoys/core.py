import os
import pickle

from .settings import STORAGE_DIR, PATHS_TO_FIND
from .utils import timeit, antifingerprint_from_paths


def run(source, target, antidecoys_filter=None, use_paths_antifp=True, verbose=True):
    filter_antifp_path = os.path.join(STORAGE_DIR, 'antifingerprint_filter.pickle')
    paths_antifp_path = os.path.join(STORAGE_DIR, 'antifingerprint_paths.pickle')
    paths_path = os.path.join(STORAGE_DIR, 'paths.pickle')
    from .pathfinders import BidirectionalPathFinder

    paths = []
    paths_antifp = None
    if os.path.exists(paths_path):
        pickled_file = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_file))
        pickled_file.close()
    if antidecoys_filter and os.path.exists(filter_antifp_path):
        pickled_file = open(filter_antifp_path, mode='rb')
        antidecoys_filter.antifingerprint = pickle.load(pickled_file)
    if use_paths_antifp and os.path.exists(paths_antifp_path):
        pickled_file = open(paths_antifp_path, mode='rb')
        paths_antifp = pickle.load(pickled_file)
        pickled_file.close()

    pathfinder = BidirectionalPathFinder(source, target, verbose=verbose, paths_antifingerprint=paths_antifp)
    for i in range(PATHS_TO_FIND):
        # find a path
        exec_time = timeit(pathfinder)
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(exec_time, i + 1))

        # reset the pathfinder for another search
        if (antidecoys_filter or use_paths_antifp) and pathfinder.path:
            # compute and save new antifingerprint
            if antidecoys_filter:
                antidecoys_filter.antifingerprint = antifingerprint_from_paths(pathfinder.path, antidecoys_filter.antifingerprint)
            if use_paths_antifp:
                paths_antifp = antifingerprint_from_paths(pathfinder.path, paths_antifp)
                pathfinder.paths_antifingerprint = paths_antifp
            pathfinder.reset()
        else:
            pathfinder = BidirectionalPathFinder(source, target, verbose=verbose, paths_antifingerprint=paths_antifp)

        # pickle the results for future use
        if antidecoys_filter:
            pickled_file = open(filter_antifp_path, mode='wb')
            pickle.dump(antidecoys_filter.antifingerprint, pickled_file)
            pickled_file.close()
        if use_paths_antifp:
            pickled_file = open(paths_antifp_path, mode='wb')
            pickle.dump(paths_antifp, pickled_file)
            pickled_file.close()
        pickled_file = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_file)
        pickled_file.close()