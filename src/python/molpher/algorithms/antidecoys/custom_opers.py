from molpher.core.operations import TreeOperation

from .utils import compute_antifp_scores

class TopScoringFilter(TreeOperation):

    def __init__(self, scores, common_bits_thrs, minimum_accept):
        super(TopScoringFilter, self).__init__()
        self.scores = scores
        self.common_bits_thrs = common_bits_thrs
        self.minimum_accept = minimum_accept

    def __call__(self):
        candidates_mask = list(self.tree.candidates_mask)
        candidates = self.tree.candidates
        for idx, mol in enumerate(candidates):
            score = self.scores[mol.smiles]
            if score > self.common_bits_thrs:
                candidates_mask[idx] = False
        accepted_count = sum(candidates_mask)
        print('Antidecoys filter identified {0} morphs below threshold.'.format(accepted_count))
        if accepted_count < self.minimum_accept:
            to_get = self.minimum_accept - accepted_count
            counter = 0
            for i in range(accepted_count, len(candidates_mask)):
                candidates_mask[i] = True
                counter += 1
                if counter >= to_get:
                    break
        print('Antidecoys filter accepted {0} morphs.'.format(sum(candidates_mask)))
        self.tree.candidates_mask = candidates_mask

class GatherAntiFPScores(TreeOperation):

    def __init__(self, out_var, anti_fp, settings):
        super(GatherAntiFPScores, self).__init__()
        self.out = out_var
        self.antitfp = anti_fp
        self.settings = settings

    def __call__(self):
        results = compute_antifp_scores(
            self.tree.candidates
            , self.antitfp
            , self.settings
        )
        self.out.update(results)