import multiprocessing
from multiprocessing.managers import BaseManager

from multiprocessing import Manager

from molpher.core.operations import TreeOperation

from .settings import *

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

class FindClosest:

    def __init__(self):
        self.closest = None

    def __call__(self, morph):
        if not self.closest:
            self.closest = morph.copy()
            return
        current_dist = self.closest.getDistToTarget()
        morph_dist = morph.getDistToTarget()
        if morph_dist < current_dist:
            self.closest = morph.copy()

class AntidecoysFilter(TreeOperation):

    def __init__(self):
        super(AntidecoysFilter, self).__init__()
        self.antifingerprint = None

    def __call__(self):
        if self.antifingerprint:
            tree = self.tree
            mask = list(tree.candidates_mask)
            candidates = [Chem.MolFromSmiles(mol.smiles) if mask[idx] else None for idx, mol in enumerate(tree.candidates)]
            candidates_fps = [Generate.Gen2DFingerprint(mol, SIG_FAC) if mol else None for mol in candidates]

            counter = 0
            for idx, fp in enumerate(candidates_fps):
                if fp:
                    candidate_bits = fp.GetNumOnBits()
                    common_fp = self.antifingerprint & fp
                    common_bits = common_fp.GetNumOnBits()
                    common_bits_perc = common_bits / candidate_bits
                    if common_bits_perc > COMMON_BITS_PERC_THRS:
                        mask[idx] = False
                        counter += 1
            print('Anti-decoys filter eliminated {0}/{1} molecules...'.format(counter, sum(mask)))
            tree.candidates_mask = mask

class AntidecoysFilterMulti(TreeOperation):

    # class ComparisonProxy(object):
    #     def __init__(self, mask):
    #         self._mask = mask
    #         self._counter = 0
    #
    #     def update_mask(self, idx, val):
    #         self._mask[idx] = val
    #         self._counter += 1
    #
    #     def get_mask(self):
    #         return self._mask
    #
    #     def get_counter(self):
    #         return self._counter
    #
    # class ComparisonManager(BaseManager):
    #     pass

    # @staticmethod
    # def get_manager():
    #     m = AntidecoysFilterMulti.ComparisonManager()
    #     m.start()
    #     return m

    @staticmethod
    def eval_morph(data):
        if data:
            mol, idx, antifingerprint, _shared = data
            fp = Generate.Gen2DFingerprint(mol, SIG_FAC)
            candidate_bits = fp.GetNumOnBits()
            common_fp = antifingerprint & fp
            common_bits = common_fp.GetNumOnBits()
            common_bits_perc = common_bits / candidate_bits
            if common_bits_perc > COMMON_BITS_PERC_THRS:
                # proxy.update_mask(idx, False)
                _shared['mask'][idx] = False
                _shared['counter'] += 1

    def __init__(self):
        super(AntidecoysFilterMulti, self).__init__()
        self.antifingerprint = None
        # self.manager = self.get_manager()
        # self.manager.register('ComparisonProxy', self.ComparisonProxy)

    def __call__(self):
        if self.antifingerprint:
            tree = self.tree
            mask = list(tree.candidates_mask)
            counter = 0

            manager = Manager()
            shared_data = {
                'mask' : mask
                , 'counter' : counter
            }
            _shared_data = manager.dict(shared_data)

            candidates = [(Chem.MolFromSmiles(mol.smiles), idx, self.antifingerprint, _shared_data) if mask[idx] else None for idx, mol in enumerate(tree.candidates)]

            pool = multiprocessing.Pool(MAX_THREADS)
            pool.map(func=self.eval_morph, iterable=candidates)
            pool.close()
            pool.join()

            print('Anti-decoys filter eliminated {0}/{1} molecules.'.format(_shared_data['counter'], sum(mask)))
            tree.candidates_mask = _shared_data['mask']