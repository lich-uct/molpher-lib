VERSION = None
BUILD_NUMBER = None

def load_stuff(path):
    with open(path, 'r') as versionfile:
        return versionfile.read().strip()

VERSION = load_stuff("VERSION.TXT")
BUILD_NUMBER = load_stuff("BUILD.TXT")

print('Version: ', VERSION)
print('Build: ', BUILD_NUMBER)