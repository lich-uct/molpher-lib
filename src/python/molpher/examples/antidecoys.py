import sys

from molpher.algorithms.antidecoys.run import run
from molpher.algorithms.antidecoys.settings import AntidecoysSettings

def main(args):
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'
    storage_dir = None
    if len(args) == 2:
        storage_dir = args[1]

    settings = AntidecoysSettings()
    settings.storage_dir = storage_dir

    run(
        cocaine
        , procaine
        , verbose=False
        , settings=settings
    )

if __name__ == "__main__":
    exit(main(sys.argv))