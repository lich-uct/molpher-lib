VERSION = None

with open('VERSION.TXT', 'r', encoding='utf-8') as versionfile:
    VERSION = versionfile.read().strip()

print(VERSION)