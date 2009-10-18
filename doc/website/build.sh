#! /bin/sh

BUILD_DIR=build
HTML_BUILD_DIR="$BUILD_DIR"/html

if [[ -z $1 ]]
then
    echo "Building full website ..."
    
    if [[ -d $HTML_BUILD_DIR ]]
    then
        echo '- Skipping creating directory "'$HTML_BUILD_DIR'": already exists.'
    else
        mkdir -p $HTML_BUILD_DIR
    fi
    
    echo '- Copying images ...'
    cp *.gif $HTML_BUILD_DIR
    cp *.png $HTML_BUILD_DIR
    
    echo '- Generating HTML ...'
    rst2html.py --stylesheet-path=dendropy.css index.rst > $HTML_BUILD_DIR/index.html
    rst2html.py --stylesheet-path=dendropy.css sumtrees.rst > $HTML_BUILD_DIR/sumtrees.html
    rst2html.py --stylesheet-path=dendropy.css tutorial.rst > $HTML_BUILD_DIR/tutorial.html
    
    echo '- Creating archive ...'
    zip -j $BUILD_DIR/dendropy-html.zip $HTML_BUILD_DIR/*
else
    rst2html.py --stylesheet-path=dendropy.css $1
fi    
