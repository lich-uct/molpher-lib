import sys

from molpher.algorithms.antidecoys.run import run
from molpher.algorithms.antidecoys.settings import AntidecoysSettings

def main(args):
    cocaine = 'CN1C2CCC1C(C(=O)OC)C(OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'
    storage_dir = None
    if len(args) == 2:
        storage_dir = args[1]
    else:
        storage_dir = 'antidecoys_data'

    settings = AntidecoysSettings(
        source=cocaine
        , target=procaine
        , storage_dir=storage_dir
        , max_threads=4
        , common_bits_max_thrs = 0.75
        , min_accepted=500
        , antidecoys_max_iters=30
        , distance_thrs=0.3
    )

    run(
        settings
        , paths_count=50
    )

if __name__ == "__main__":
    exit(main(sys.argv))