import warnings

import molpher.swig_wrappers.core as wrappers

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

    @property
    def smiles(self):
        return self.getSMILES()

    @smiles.setter
    def smiles(self, val):
        self.setSMILES(val)

    @property
    def tree(self):
        return self.getTree()