import os

import multiprocessing


class Settings:

    def __init__(
            self
            , source
            , target
            , storage_dir
            , max_threads = None
            , max_iters = None
            , tree_params = None
            , verbose = False
    ):
        self.verbose = verbose

        self.source = source
        self.target = target

        self.storage_dir = storage_dir
        if not os.path.exists(self.storage_dir):
            os.mkdir(self.storage_dir)

        self.max_threads = max_threads if max_threads not in (None, False) else multiprocessing.cpu_count()
        self.max_iters = max_iters if max_iters not in (None, False) else 100
        self.tree_params = tree_params