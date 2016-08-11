import os

import pickle

from .pathfinder import BidirectionalPathFinder


def run(
        source
        , target
        , storage_dir
        , verbose=False
        , paths_to_find=1
):
    paths_path = os.path.join(storage_dir, 'paths.pickle')
    if not os.path.exists(storage_dir):
        os.mkdir(storage_dir)

    paths = []
    if os.path.exists(paths_path):
        pickled_paths = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_paths))
        pickled_paths.close()
    for i in range(paths_to_find):
        # find a path
        pathfinder = BidirectionalPathFinder(source, target, verbose=verbose)
        pathfinder()
        paths.append(pathfinder.path)

        # pickle the results for future use
        pickled_paths = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_paths)
        pickled_paths.close()