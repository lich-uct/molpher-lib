"""
Custom operations and callbacks shared between the algorithms.

"""

class FindClosest:
    """
    Callback that will search for a molecule closest to the target and
    saves it as a copy to the self.closest

    .. seealso:: :class:`~molpher.core.operations.TraverseOper`, :meth:`~molpher.core.ExplorationTree.ExplorationTree.traverse`

    """

    def __init__(self):
        self.closest = None
        """:class:`~molpher.core.MolpherMol.MolpherMol` instance (copy of the closest molecule)
        after the callback was executed, initialized to `None`"""

    def __call__(self, morph):
        if not self.closest:
            self.closest = morph.copy()
            return
        current_dist = self.closest.getDistToTarget()
        morph_dist = morph.getDistToTarget()
        if morph_dist < current_dist:
            self.closest = morph.copy()