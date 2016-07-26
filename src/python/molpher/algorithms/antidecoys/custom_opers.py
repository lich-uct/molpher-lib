import multiprocessing

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
                    assert common_bits_perc <= 1.0
                    if common_bits_perc > COMMON_BITS_PERC_THRS:
                        mask[idx] = False
                        counter += 1
            print('Anti-decoys filter eliminated {0}/{1} molecules...'.format(counter, sum(mask)))
            tree.candidates_mask = mask

class AntidecoysFilterMulti(TreeOperation):

    @staticmethod
    def eval_morph(data):
        if data:
            mol, idx, antifingerprint, shared_list, shared_counter = data
            fp = Generate.Gen2DFingerprint(mol, SIG_FAC)
            candidate_bits = fp.GetNumOnBits()
            common_fp = antifingerprint & fp
            common_bits = common_fp.GetNumOnBits()
            common_bits_perc = common_bits / candidate_bits
            assert common_bits_perc <= 1.0
            if common_bits_perc > COMMON_BITS_PERC_THRS:
                # proxy.update_mask(idx, False)
                shared_list[idx] = False
                shared_counter.value += 1

    def __init__(self):
        super(AntidecoysFilterMulti, self).__init__()
        self.antifingerprint = None

    def __call__(self):
        if self.antifingerprint:
            tree = self.tree
            mask = list(tree.candidates_mask)

            manager = Manager()
            # shared_data = {
            #     'counter' : counter
            # }
            _shared_list = manager.list(mask)
            _shared_counter = manager.Value('I', value=0)

            candidates = [(Chem.MolFromSmiles(mol.smiles), idx, self.antifingerprint, _shared_list, _shared_counter) if mask[idx] else None for idx, mol in enumerate(tree.candidates)]

            pool = multiprocessing.Pool(MAX_THREADS)
            pool.map(func=self.eval_morph, iterable=candidates)
            pool.close()
            pool.join()

            print('Anti-decoys filter eliminated {0}/{1} molecules.'.format(_shared_counter.value, sum(mask)))
            tree.candidates_mask = _shared_list
