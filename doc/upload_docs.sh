#!/usr/bin/env bash

set -e

COMMIT_ID=`git log -1 --pretty=short --abbrev-commit`
export VERSION=`cat ../VERSION.TXT`
MSG="Adding docs to gh-pages for $COMMIT_ID (v$VERSION)"

BASE_DIR="`pwd`"
HTML_DIR="$BASE_DIR/build/html/"

TMPREPO=/tmp/docs/$USER/
rm -rf $TMPREPO
mkdir -p -m 0755 $TMPREPO
git clone `git config --get remote.origin.url` $TMPREPO

cd $TMPREPO
git checkout gh-pages
cp $BASE_DIR/homepage/* $TMPREPO
python render.py
rm -f *.py
rm -f *.template
mkdir -p "v$VERSION/"
cp -r $HTML_DIR/* "v$VERSION/"
touch .nojekyll

git add -A
git commit -m "$MSG" && git push origin gh-pages