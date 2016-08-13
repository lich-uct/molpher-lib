import os

import pickle

from .pathfinder import BidirectionalPathFinder


def run(
        settings
        , paths_to_find=1
):
    paths_path = os.path.join(settings.storage_dir, 'paths.pickle')

    paths = []
    if os.path.exists(paths_path):
        pickled_paths = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_paths))
        pickled_paths.close()
    for i in range(paths_to_find):
        # find a path
        pathfinder = BidirectionalPathFinder(settings)
        pathfinder()
        paths.append(pathfinder.path)

        # pickle the results for future use
        pickled_paths = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_paths)
        pickled_paths.close()

    return paths