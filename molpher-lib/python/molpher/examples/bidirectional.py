import sys

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols

ACCEPT_MAX = 10

class MyFilterMorphs(TreeOperation):

    def __init__(self):
        super(MyFilterMorphs, self).__init__()

    def __call__(self):
        tree = self.getTree()
        mask = [False for x in tree.candidates_mask]
        for i in range(ACCEPT_MAX):
            mask[i] = True
        tree.candidates_mask = mask

    def getTree(self):
        tree = super(MyFilterMorphs, self).getTree()
        if tree:
            tree.__class__ = ETree # 'cast' the wrapped class to the 'pretty' Python proxy class
        return tree

class AdjustDistances(TreeOperation):

        def __init__(self):
            super(AdjustDistances, self).__init__()
            self.other = None
            self.path_found = False
            self.minimum_distance = sys.float_info.max
            self.closest_pair = None
            self.max_closest_to_adjust = ACCEPT_MAX

        def __call__(self):
            super(AdjustDistances, self).__call__()
            if self.tree.params['source'] != self.other.params['source'] and self.other:
                other_morphs = sorted([x for x in self.other.candidates], key=lambda x : x.getDistToTarget())[:self.max_closest_to_adjust]
                this_morphs = sorted([x for x in self.tree.candidates], key=lambda x : x.getDistToTarget())[:self.max_closest_to_adjust]
                print('Closest to target (this): {0}'.format(this_morphs[0].getDistToTarget()))
                print('Closest to target (other): {0}'.format(other_morphs[0].getDistToTarget()))
                assert other_morphs and this_morphs

                other_fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x.getSMILES())) for x in other_morphs]
                this_fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x.getSMILES())) for x in this_morphs]

                morph_other_smiles_min = None
                for idx_this, morph_this in enumerate(this_morphs):
                    min_dist = sys.float_info.max
                    morph_this_smiles = morph_this.getSMILES()
                    morph_other_smiles = None
                    min_idx_other = 0
                    for idx_other, morph_other in enumerate(other_morphs):
                        morph_other_smiles = morph_other.getSMILES()
                        dist = DataStructs.FingerprintSimilarity(this_fps[idx_this], other_fps[idx_other])
                        if dist < min_dist:
                            min_dist = dist
                            min_idx_other = idx_other
                            morph_other_smiles_min = morph_other_smiles

                            if morph_this_smiles == morph_other_smiles or dist == 0:
                                print('Path Found:')
                                print('Two equal candidates are: {0} and {1}'.format(morph_this_smiles, morph_other_smiles))
                                self.path_found = True

                    target_dist = morph_this.getDistToTarget()
                    if min_dist < target_dist:
                        morph_this.setDistToTarget(min_dist)
                        other_morphs[min_idx_other].setDistToTarget(min_dist)
                    if min_dist < self.minimum_distance:
                        self.minimum_distance = min_dist
                        self.closest_pair = (morph_this_smiles, morph_other_smiles_min)

            else:
                raise AttributeError('Error: Same tree specified as other or no tree specified at all.')


        def setOther(self, other):
            self.other = other

class BidirectionalPathFinder:

    def __init__(self, source, target):
        self.source_target = ETree(source=source, target=target)
        self.target_source = ETree(source=target, target=source)
        self.distance_adjust = AdjustDistances()
        self.distance_adjust.setOther(self.target_source)
        self.ITERATION = [
            GenerateMorphsOper()
            , self.distance_adjust
            , SortMorphsOper()
            , MyFilterMorphs()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

    def __call__(self):
        counter = 0
        while True:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                if oper.__class__ != AdjustDistances:
                    self.target_source.runOperation(oper)
                self.source_target.runOperation(oper)

            print('Closest pair:', self.distance_adjust.closest_pair)
            print('Distance:', self.distance_adjust.minimum_distance)

            if self.distance_adjust.path_found:
                break

cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

pathfinder = BidirectionalPathFinder(cocaine, procaine)
pathfinder()
