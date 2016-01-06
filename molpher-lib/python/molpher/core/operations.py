"""
This module just contains references to
the operations defined
in `molpher.swig_wrappers.core`.

"""

import molpher
from molpher.core import ExplorationTree
from molpher.swig_wrappers.core import \
    GenerateMorphsOper \
    , FilterMorphsOper \
    , ExtendTreeOper \
    , PruneTreeOper \
    , SortMorphsOper \
    , FindLeavesOper \
    , TraverseOper \
    , TraverseCallback

class TreeOperation(molpher.swig_wrappers.core.TreeOperation):

    def __init__(self):
        super(TreeOperation, self).__init__()
        self.tree = self.getTree()

    def __call__(self):
        self.tree = self.getTree()

    def getTree(self):
        tree = super(TreeOperation, self).getTree()
        if tree:
            tree.__class__ = ExplorationTree # 'cast' the wrapped class to the 'pretty' Python proxy class
        return tree