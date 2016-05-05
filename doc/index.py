import os
import shutil

from jinja2 import Template

# The short X.Y version.
current_version = 'v' + os.environ['VERSION']
versions = [ x for x in os.listdir(".") if x.startswith("v")]

if current_version in versions:
    print("WARNING: Docs for this version already built. Rewriting...")
    shutil.rmtree(current_version)
else:
    versions.append(current_version)

INDEX_TEMPLATE_PATH = "index.html.template"
INDEX_PATH = "index.html"

versions.sort()

template = None
with open(INDEX_TEMPLATE_PATH, "r", encoding="utf-8") as template_file:
    template = Template(template_file.read())

with open(INDEX_PATH, "w", encoding='utf-8') as index_file:
    index_file.write(template.render(
        versions=versions
    ))