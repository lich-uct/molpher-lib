from molpher.algorithms.classic.run import run


def main():
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    path = run(cocaine, procaine)
    print('Path found: {0}'.format(path))

if __name__ == "__main__":
    exit(main())