from rdkit import Chem
from molpher.core import MolpherMol
from molpher.core.morphing import Molpher
from molpher.core.morphing.operators import *

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
    ]
    , attempts = 100 # create at most 100 molecules
    , collectors = [collect_sensible]
)

# execute morphing and show created molecules
molpher()
as_mol_grid(sensible_morphs.values()) # draw generated structures in a grid
