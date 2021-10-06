from rdkit import Chem
from molpher.core import MolpherMol, MolpherAtom, ExplorationTree as ETree
from molpher.core.morphing.operators import MorphingOperator
from molpher.random_numbers import get_random_number

class AddFragment(MorphingOperator):

    def __init__(self, fragment, open_atoms_frag, oper_name):
        super(AddFragment, self).__init__()
        self._name = oper_name
        self._fragment = fragment
        self._open_atoms_frag = open_atoms_frag
        self._orig_rdkit = None
        self._open_atoms = []

    def setOriginal(self, mol):
        super(AddFragment, self).setOriginal(mol)
        if self.original:
            self._orig_rdkit = self.original.asRDMol()
            self._open_atoms = []

            for atm_rdkit, atm_molpher in zip(self._orig_rdkit.GetAtoms(), self.original.atoms):
                free_bonds = atm_rdkit.GetImplicitValence()
                if free_bonds >= 1 and not (MolpherAtom.NO_ADDITION & atm_molpher.locking_mask):
                    self._open_atoms.append(atm_rdkit.GetIdx())

    def morph(self):
        combo_mol = Chem.EditableMol(Chem.CombineMols(
            self._orig_rdkit
            , self._fragment
        ))
        atom_orig = self._open_atoms[get_random_number(0, len(self._open_atoms)-1)]
        atom_frag = len(self.original.atoms) + self._open_atoms_frag[get_random_number(0, len(self._open_atoms_frag)-1)]
        combo_mol.AddBond(atom_orig, atom_frag, order=Chem.rdchem.BondType.SINGLE)
        combo_mol = combo_mol.GetMol()
        Chem.SanitizeMol(combo_mol)

        ret = MolpherMol(other=combo_mol)
        for atm_ret, atm_orig in zip(ret.atoms, self.original.atoms):
            atm_ret.locking_mask = atm_orig.locking_mask

        return ret

    def getName(self):
        return self._name

def main(captopril=None):
    mol = MolpherMol("CC=O")
    frag = Chem.MolFromSmiles('c1ccccc1')
    oper = AddFragment(frag, [1], "Add Benzyl")
    oper.setOriginal(mol)
    morph = oper.morph()
    print(morph.smiles)

    if not captopril:
        captopril = MolpherMol("src/python/molpher/examples/captopril.sdf")
    tree = ETree.create(source=captopril)
    tree.morphing_operators = tree.morphing_operators + (oper,)
    print(tree.morphing_operators)

    tree.generateMorphs()
    print([x.smiles for x in tree.candidates])

if __name__ == '__main__':
   main()