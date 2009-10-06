#! /bin/sh

if [[ -d html ]]
then
    echo '- Skipping creating directory "html": already exists.'
else
    mkdir html
fi

echo '- Copying images ...'
cp *.gif html

echo '- Generating HTML ...'
rst2html.py sumtrees.rst > html/sumtrees.html

echo '- Creating archive ...'
zip -j dendropy-html.zip html/*

