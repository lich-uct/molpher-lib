import time

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds


def find_path(tree, end_mol):
    path = []
    current = tree.fetchMol(end_mol)
    path.append(current.getSMILES())
    while current != '':
        current = current.getParentSMILES()
        if current:
            current = tree.fetchMol(current)
            path.append(current.getSMILES())
    path.reverse()
    return path