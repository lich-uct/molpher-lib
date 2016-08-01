import multiprocessing

from multiprocessing import Manager

from molpher.core.operations import TreeOperation

from .settings import SIG_FAC, COMMON_BITS_PERC_THRS, MAX_THREADS
from .utils import compute_antifp_scores

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

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
            shared_mask = manager.list(mask)
            shared_counter = manager.Value('I', value=0)

            candidates_data = [
                (Chem.MolFromSmiles(mol.smiles), idx, self.antifingerprint, shared_mask, shared_counter)
                if mask[idx]
                else None
                for idx, mol in enumerate(tree.candidates)
                ]

            pool = multiprocessing.Pool(MAX_THREADS)
            pool.map(func=self.eval_morph, iterable=candidates_data)
            pool.close()
            pool.join()

            tree.candidates_mask = shared_mask
            print('Anti-decoys filter eliminated {0}/{1} molecules. Remaining: {2}.'.format(shared_counter.value, sum(mask), sum(tree.candidates_mask)))

class TopXFilter(TreeOperation):

    def __init__(self, top_max_count):
        super(TopXFilter, self).__init__()
        self.top_max_count = top_max_count

    def __call__(self):
        candidates_mask = list(self.tree.candidates_mask)
        if len(candidates_mask) > self.top_max_count:
            for i in range(self.top_max_count):
                candidates_mask[i] = True
            for i in range(self.top_max_count, len(candidates_mask)):
                candidates_mask[i] = False
        self.tree.candidates_mask = candidates_mask

class GatherAntiFPScores(TreeOperation):

    def __init__(self, out_var, anti_fp):
        super(GatherAntiFPScores, self).__init__()
        self.out = out_var
        self.antitfp = anti_fp

    def __call__(self):
        results = compute_antifp_scores(self.tree.candidates, self.antitfp)
        self.out.update(results)