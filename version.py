import os

VERSION = None
BUILD_NUMBER = None

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def load_versionfile(path):
    with open(path, 'r') as versionfile:
        return versionfile.read().strip()

VERSION = load_versionfile(os.path.join(BASE_DIR, "VERSION.TXT"))
BUILD_NUMBER = load_versionfile(os.path.join(BASE_DIR, "BUILD.TXT"))

print('Version: ', VERSION)
print('Build: ', BUILD_NUMBER)