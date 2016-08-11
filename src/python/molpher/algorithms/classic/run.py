from molpher.algorithms.classic.pathfinder import ClassicPathFinder


def run(source, target):

    pathfinder = ClassicPathFinder(source, target)
    pathfinder()
    return pathfinder.path