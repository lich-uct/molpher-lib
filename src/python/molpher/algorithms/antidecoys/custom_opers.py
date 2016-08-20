"""
Custom operations and callbacks used in
`AntidecoysPathFinder`.

"""

from molpher.core.operations import TreeOperation
from molpher.core.operations.callbacks import SortMorphsCallback

from .utils import compute_antifp_scores

class GatherAntiFPScores(TreeOperation):
    """
    :param out: a dictionary to save the results to
    :param anti_fp: the anti-fingerprint
    :param settings: `AntidecoysSettings`

    A `tree operation` which computes the anti-fingerprint scores using the
    given settings and anti-fingerprint. For every
    molecule in `candidates` the similarity
    to the anti-fingerprint is computed and
    saved as a score into the supplied dictionary.

    ..  seealso:: `compute_antifp_scores`

    """

    def __init__(self, out, anti_fp, settings):
        super(GatherAntiFPScores, self).__init__()
        self.out = out
        """output dictionary"""
        self.antitfp = anti_fp
        """the anti-fingerprint"""
        self.settings = settings
        """`AntidecoysSettings`"""

    def __call__(self):
        results = compute_antifp_scores(
            self.tree.candidates
            , self.antitfp
            , self.settings
        )
        self.out.update(results)

class AntiFpSortCallback(SortMorphsCallback):
    """
    :param antifp_scores: SMILES to anti-fingerprint similarity map (for each molecule in `candidates`)
    :type antifp_scores: `dict`

    This callback is used to sort `candidates` with regard
    to their similarity to the anti-fingerprint. Molecules
    less similar to the anti-fingerprint are ranked higher
    than more similar ones.

    """

    def __init__(self, antifp_scores):
        super(AntiFpSortCallback, self).__init__()
        self.minimum_common_bits_perc = 1.0
        """minimum percentage of bits that one
        `candidate morph <candidate morphs>` has in common with the anti-fingerprint"""
        self.maximum_common_bits_perc = 0.0
        """maximum percentage of bits that one
        `candidate morph <candidate morphs>` has in common with the anti-fingerprint"""
        self.antifp_scores = antifp_scores
        """SMILES to anti-fingerprint score map (for each molecule in `candidates`)"""

    def __call__(self, a, b):
        """
        Overrides the abstract method in the base class.
        Compares two morphs according to the percentage
        of bits that they have in common with the
        anti-fingerprint.

        :param a: morph on the left side of the comparison
        :param b: morph on the right side of the comparison
        :return: `True` if morph on the right side is less similar to the anti-fingerprint than the one on the left, `False` otherwise
        """

        perc_a = self.antifp_scores[a.getSMILES()]
        perc_b = self.antifp_scores[b.getSMILES()]
        minimum = min(perc_a, perc_b, self.minimum_common_bits_perc)
        maximum = max(perc_a, perc_b, self.minimum_common_bits_perc)
        if minimum < self.minimum_common_bits_perc:
            self.minimum_common_bits_perc = minimum
        if maximum > self.maximum_common_bits_perc:
            self.maximum_common_bits_perc = maximum
        return perc_a < perc_b

class AntidecoysFilter(TreeOperation):
    """
    :param scores: a `dict` mapping a SMILES string of a molecule in `candidates` to the corresponding score value
    :param common_bits_thrs: maximum percentage of bits in common with the anti-fingerprint allowed in the
    :param minimum_accept: always accept at least this many top scoring molecules

    A filtering `tree operation` which filters
    the morphs in `candidates` based
    on :samp:`common_bits_thrs`. All morphs with
    anti-fingerprint score bigger than
    :samp:`common_bits_thrs` are rejected. However,
    at least :samp:`minimum_accept` best morphs
    survives.

    """

    def __init__(self, scores, common_bits_thrs, minimum_accept):
        super(AntidecoysFilter, self).__init__()
        self.scores = scores
        """a `dict` mapping a SMILES string of a molecule in `candidates` to the corresponding score value"""
        self.common_bits_thrs = common_bits_thrs
        """maximum percentage of bits in common with the anti-fingerprint allowed in the"""
        self.minimum_accept = minimum_accept
        """always accept at least this many top scoring molecules"""

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