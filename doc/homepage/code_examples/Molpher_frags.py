from rdkit import Chem
from molpher.core import MolpherMol, MolpherAtom
from molpher.core.morphing import Molpher
from molpher.core.morphing.operators import *
from molpher.random import get_random_number

class AddFragment(MorphingOperator):
    """
    Attaches a given molecule fragment to an atom in the molecule.
    """

    def __init__(self, fragment, open_atoms_frag, oper_name):
        super(AddFragment, self).__init__()
        self._name = oper_name # name of the operator
        self._fragment = fragment # fragment as RDKit Mol
        self._open_atoms_frag = open_atoms_frag # possible attachment positions on the fragment
        self._orig_rdkit = None # original molecule as RDKit Mol
        self._open_atoms = [] # possible attachment positions on the original molecule

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

# define a collector -> a callback function that processes morphs as they are generated
strange_patterns = Chem.MolFromSmarts('[S,O,N][F,Cl,Br,I]')
sensible_morphs = dict()
def collect_sensible(morph, operator):
    """
    simple collector, accepts morphs without some weird structural patterns
    """

    rd_morph = morph.asRDMol()
    if not rd_morph.HasSubstructMatch(strange_patterns):
        sensible_morphs[morph.smiles] = morph
        morph.parent_operator = operator.getName()

# create some AddFragment operators
fragments = ['c1ccccc1', 'C(=O)O']
add_frags = []
for frag in fragments:
    add_frag = AddFragment(Chem.MolFromSmiles(frag), [0], "Add " + frag)
    add_frags.append(add_frag)

# load a molecule from SDF and generate some derived molecules with given morphing operators 
mol = MolpherMol("captopril.sdf")
molpher = Molpher(
    mol
    , [ # list of morphing operators to use
        AddAtom()
        , RemoveAtom()
        , MutateAtom()
        , AddBond()
        , RemoveBond()
        , ContractBond()
        , InterlayAtom()
        , RerouteBond()
    ] + add_frags # add our custom operators, too
    , attempts = 100 # create at most 100 molecules
    , collectors = [collect_sensible]
)

# execute morphing and show created molecules
molpher()
as_mol_grid(sensible_morphs.values()) # draw generated structures in a grid
