from molpher.algorithms.classic.run import run
from molpher.algorithms.settings import Settings

# our source and target molecules
cocaine = 'CN1C2CCC1C(C(=O)OC)C(OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

# directory where the path will be stored (as a pickled list)
storage_dir = 'data'

# initialize the exploration settings
settings = Settings(
    source=cocaine
    , target=procaine
    , storage_dir=storage_dir
    , max_threads=4
)

run(settings)
