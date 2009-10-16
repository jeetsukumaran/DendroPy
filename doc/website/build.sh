#! /bin/sh

BUILD_DIR=html

if [[ -d $BUILD_DIR ]]
then
    echo '- Skipping creating directory "'$BUILD_DIR'": already exists.'
else
    mkdir $BUILD_DIR
fi

echo '- Copying images ...'
cp *.gif $BUILD_DIR

echo '- Generating HTML ...'
rst2html.py --stylesheet-path=dendropy.css index.rst > $BUILD_DIR/index.html
rst2html.py --stylesheet-path=dendropy.css sumtrees.rst > $BUILD_DIR/sumtrees.html
rst2html.py --stylesheet-path=dendropy.css cookbook.rst > $BUILD_DIR/cookbook.html

echo '- Creating archive ...'
zip -j dendropy-html.zip $BUILD_DIR/*

