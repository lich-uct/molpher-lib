import os
from shutil import copyfile

import subprocess
from jinja2 import Template
import version

BUILD_ALL_PYTHON = False
TARGETS = ["tbb", "molpher-lib"]
PYTHON_VERSIONS = ['2.6', '2.7',  '3.3',  '3.4', '3.5'] if BUILD_ALL_PYTHON else ['3.5', '2.7']
VERSION = version.VERSION
BUILD_NUMBER = version.BUILD_NUMBER
LICENSE_FILE = "LICENSE.txt"
TBB_LICENSE_FILE = "deps/tbb/COPYING"
TBB_VERSION = "4.2"

for target in TARGETS:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(BASE_DIR)
    for python_version in PYTHON_VERSIONS:
        PACKAGE_DIR = os.path.join(BASE_DIR, "conda/{0}/{0}/".format(target))
        METAFILE_TEMPLATE_PATH = os.path.join(PACKAGE_DIR, "../meta.yaml.template")
        METAFILE_PATH = os.path.join(PACKAGE_DIR, "meta.yaml")

        template = None
        with open(METAFILE_TEMPLATE_PATH, "r", encoding="utf-8") as metafile_template:
            template = Template(metafile_template.read())

        with open(METAFILE_PATH, "w", encoding='utf-8') as metafile:
            metafile.write(template.render(
                target=target
                , version=VERSION + '_r' + BUILD_NUMBER if BUILD_NUMBER != '0' else VERSION
                , build_number=BUILD_NUMBER
                , license_file=LICENSE_FILE
                , build_string="py{0}_{1}".format(python_version.replace('.', ''), BUILD_NUMBER)
                , python_spec="python {0}*".format(python_version)
                , tbb_version=TBB_VERSION
                , tbb_spec = "tbb {0}*".format(TBB_VERSION)
            ))

        os.chdir(os.path.join(PACKAGE_DIR, "../"))
        os.environ['BASE_DIR'] = BASE_DIR
        if target == 'tbb':
            # do not set Python version for tbb
            copyfile(os.path.join(BASE_DIR, TBB_LICENSE_FILE), os.path.join(PACKAGE_DIR, LICENSE_FILE))
            ret = subprocess.call("conda build {0}".format(target), env=os.environ, shell=True)
        else:
            copyfile(os.path.join(BASE_DIR, LICENSE_FILE), os.path.join(PACKAGE_DIR, LICENSE_FILE))
            ret = subprocess.call("conda build {0} --python {1}".format(target, python_version), env=os.environ, shell=True)
        if ret != 0:
            print("An error has occured during conda build...")
            exit(1)

        # remove the generated and copied files and go back to BASE_DIR
        os.remove(METAFILE_PATH)
        os.remove(os.path.join(PACKAGE_DIR, LICENSE_FILE))
        os.chdir(BASE_DIR)

        # make sure that the tbb package is only built once
        if target == 'tbb':
            break
