from molpher.algorithms.antidecoys.utils import run

def main():
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'
    run(cocaine, procaine)

if __name__ == "__main__":
    exit(main())