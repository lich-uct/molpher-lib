import warnings

import molpher.swig_wrappers.core as wrappers
import molpher.core.ExplorationTree


class MolpherMol(wrappers.MolpherMol):

    def __init__(self, smiles=None, other=None):
        if smiles and other:
            warnings.warn(
                "Both SMILES string and another MolpherMol instance specified. Using another instance to initialize..."
                , RuntimeWarning
            )
        if other and isinstance(other, wrappers.MolpherMol):
            super(MolpherMol, self).__init__(other)
        elif smiles and type(smiles) == str:
            super(MolpherMol, self).__init__(smiles)
        elif not smiles and not other:
            super(MolpherMol, self).__init__()
        else:
            raise AttributeError('Invalid argumetns supplied to the constructor.')

    def copy(self):
        copy = super(MolpherMol, self).copy()
        copy.__class__ = MolpherMol
        return copy

    @property
    def tree(self):
        ret = self.getTree()
        if ret:
            ret.__class__ = molpher.core.ExplorationTree
        return ret

    @property
    def smiles(self):
        return self.getSMILES()

    @smiles.setter
    def smiles(self, val):
        self.setSMILES(val)

    @property
    def parent_smiles(self):
        return self.getParentSMILES()

    @parent_smiles.setter
    def parent_smiles(self, val):
        self.setParentSMILES(val)

    @property
    def parent_operator(self):
        return wrappers.ChemOperShortDesc(self.getParentOper())

    @property
    def dist_to_target(self):
        return self.getDistToTarget()

    @dist_to_target.setter
    def dist_to_target(self, dist):
        self.setDistToTarget(dist)

    @property
    def sascore(self):
        return self.getSAScore()

    @sascore.setter
    def sascore(self, val):
        self.setSAScore(val)

    @property
    def historic_descendents(self):
        return self.getHistoricDescendants()

    @historic_descendents.setter
    def historic_descendents(self, val):
        self.setHistoricDescendants(val)